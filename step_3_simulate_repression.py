import pandas as pd
import numpy as np
import os
from scipy.optimize import curve_fit
from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d
from sklearn.metrics import mean_absolute_error
from plotnine import (
    ggplot,
    aes,
    geom_point,
    geom_line,
    labs,
    coord_cartesian,
)
from utils import (
    load_nfc_background,
    load_fcs_data,
    theme_castuner,
    color_bfp,
    color_mcherry,
    plasmid_map
)


# Define ODE model for repression
def ode_model_kd(t, y, t_up, K, n, alpha):
    R, Y = y
    # R approaches 1 with half-life t_up
    dR_dt = (np.log(2) / t_up) * (1 - R)
    # Y is repressed by R
    dY_dt = (K ** n / (K ** n + R ** n)) - alpha * Y
    return [dR_dt, dY_dt]


def main():
    OUT_PATH = "plots"
    PARAM_PATH = "parameters"
    os.makedirs(OUT_PATH, exist_ok=True)
    os.makedirs(PARAM_PATH, exist_ok=True)

    mBFP_neg, mmCherry_neg = load_nfc_background()
    medianexp = load_fcs_data("fcs_files/time-course_data", mBFP_neg, mmCherry_neg)

    parts = medianexp['filename'].str.split('_', expand=True)
    medianexp['plasmid'] = parts[2]
    medianexp['exp'] = parts[3]
    medianexp['time'] = pd.to_numeric(parts[5])

    KD = medianexp[medianexp['exp'] == 'KD'].copy()

    # Calculate mCherry fold change
    t0 = KD[KD['time'] == 0].groupby('plasmid')['PE-A'].mean().reset_index(name='mmcherry')
    KD = pd.merge(KD, t0, on='plasmid', how='left')
    KD['fc.cherry'] = KD['PE-A'] / KD['mmcherry']

    # Min-max scaling of bfp
    mean_final = KD[KD['time'] > 10].groupby('plasmid')['BV421-A'].mean().reset_index(name='mean_final')
    mean_init = KD[KD['time'] == 0].groupby('plasmid')['BV421-A'].mean().reset_index(name='mean_init')
    KD = pd.merge(KD, mean_final, on='plasmid', how='left')
    KD = pd.merge(KD, mean_init, on='plasmid', how='left')
    KD['norm.bfp'] = (KD['BV421-A'] - KD['mean_init']) / (KD['mean_final'] - KD['mean_init'])

    # --- Load all parameters ---
    t_down = pd.read_csv(os.path.join(PARAM_PATH, 'half_times_downregulation.csv')).rename(
        columns={'halftime': 't_down'})
    t_up = pd.read_csv(os.path.join(PARAM_PATH, 'half_times_upregulation.csv')).rename(columns={'halftime': 't_up'})
    hill_pars = pd.read_csv(os.path.join(PARAM_PATH, 'Hill_parameters.csv'))
    alpha_df = pd.read_csv(os.path.join(PARAM_PATH, 'alphamcherry.csv'))
    alpha = alpha_df['alpha'].iloc[0]

    all_par = pd.merge(t_down, t_up, on='plasmid')
    all_par = pd.merge(all_par, hill_pars, on='plasmid')
    all_par['alpha'] = alpha

    # --- Simulate ODEs and find delays ---
    plasmids_to_run = ['SP430ABA', 'SP428', 'SP427', 'SP430', 'SP411']
    delays_kd = []

    for data_plasmid in plasmids_to_run:
        param_plasmid = plasmid_map[data_plasmid]

        data = KD[KD['plasmid'] == data_plasmid]
        if data.empty:
            print(f"No data for {data_plasmid}, skipping.")
            continue

        try:
            sel_par = all_par[all_par['plasmid'] == param_plasmid].iloc[0]
        except IndexError:
            print(f"No parameters for {param_plasmid}, skipping.")
            continue

        # Initial conditions
        # R starts at 0 (min normalized bfp)
        # Y starts at 1/alpha (max mCherry)
        y0_ode = [0, 1 / sel_par['alpha']]

        t_span = [0, 150]
        t_eval = np.linspace(0, 150, 300)

        exp_data = data.groupby('time')['fc.cherry'].mean().reset_index()

        # Test delays
        delay_times = np.arange(0, 25.5, 0.5)
        mae_results = []

        for t_delay in delay_times:
            sol = solve_ivp(
                ode_model_kd,
                t_span,
                y0_ode,
                t_eval=t_eval,
                args=(sel_par['t_up'], sel_par['K'], sel_par['n'], sel_par['alpha']),
                method='LSODA'
            )

            sim_time = sol.t + t_delay
            sim_Y_fc = sol.y[1] * sel_par['alpha']

            f_interp = interp1d(sim_time, sim_Y_fc, bounds_error=False, fill_value="extrapolate")
            sim_at_exp_time = f_interp(exp_data['time'])

            valid_idx = ~np.isnan(sim_at_exp_time) & ~np.isnan(exp_data['fc.cherry'])
            if np.sum(valid_idx) > 1:
                # R script calculates MAE as sum(abs(res)) / (N-1)
                N = np.sum(valid_idx)
                mae = np.sum(np.abs(exp_data['fc.cherry'][valid_idx] - sim_at_exp_time[valid_idx])) / (N - 1)
                mae_results.append({'t': t_delay, 'MAE': mae})

        mae_df = pd.DataFrame(mae_results)
        best_delay = mae_df.loc[mae_df['MAE'].idxmin()]
        delays_kd.append({'plasmid': param_plasmid, 'd_rev': best_delay['t']})

        # --- Plot MAE ---
        p_mae = (
                ggplot(mae_df, aes('t', 'MAE'))
                + geom_point(size=0.1, alpha=0.4, color="black")
                + geom_point(
            data=best_delay.to_frame().T,
            aes(x='t', y='MAE'),
            size=0.8,
            alpha=1,
            color=color_mcherry
        )
                + geom_line(aes(y=best_delay['MAE']), linetype='dashed')
                + coord_cartesian(xlim=[0, 25], ylim=[0, 0.3])
                + labs(y="MAE", x="Delay (hours)")
                + theme_castuner()
        )
        p_mae.save(
            os.path.join(OUT_PATH, f"MAE_KD_{param_plasmid}_mcherry.pdf"),
            width=1.5 * 1.618,
            height=1.5,
            units="in"
        )

        # --- Plot final simulation ---
        sol_no_delay = solve_ivp(
            ode_model_kd, t_span, y0_ode, t_eval=t_eval,
            args=(sel_par['t_up'], sel_par['K'], sel_par['n'], sel_par['alpha']),
            method='LSODA'
        )
        sim_no_delay_df = pd.DataFrame({
            'time': sol_no_delay.t,
            'R': sol_no_delay.y[0],
            'Y': sol_no_delay.y[1] * sel_par['alpha']
        })

        sim_delay_df = pd.DataFrame({
            'time': sol_no_delay.t + best_delay['t'],
            'R': sol_no_delay.y[0],
            'Y': sol_no_delay.y[1] * sel_par['alpha']
        })

        p_mcherry = (
                ggplot(data, aes('time', 'fc.cherry'))
                + geom_point(size=0.8, alpha=0.4, color=color_mcherry)
                + geom_line(data=sim_no_delay_df, aes(x='time', y='Y'), color='black', linetype='solid')
                + geom_line(data=sim_delay_df, aes(x='time', y='Y'), color='black', linetype='dashed')
                + coord_cartesian(xlim=[0, 150], ylim=[0, 1.3])
                + labs(x="Time (hours)", y="mCherry")
                + theme_castuner()
        )
        p_mcherry.save(
            os.path.join(OUT_PATH, f"KD_ODE_mCherry_{param_plasmid}.pdf"),
            width=1.5 * 1.618,
            height=1.5,
            units="in"
        )

        p_bfp = (
                ggplot(data, aes('time', 'norm.bfp'))
                + geom_point(size=0.8, alpha=0.4, color=color_bfp)
                + geom_line(data=sim_no_delay_df, aes(x='time', y='R'), color='black', linetype='solid')
                + coord_cartesian(xlim=[0, 150], ylim=[0, 1.3])
                + labs(x="Time (hours)", y="tagBFP")
                + theme_castuner()
        )
        p_bfp.save(
            os.path.join(OUT_PATH, f"KD_ODE_tagBFP_{param_plasmid}.pdf"),
            width=1.5 * 1.618,
            height=1.5,
            units="in"
        )

    delays_df = pd.DataFrame(delays_kd)
    delays_df.to_csv(
        os.path.join(PARAM_PATH, "delays_repression.csv"),
        index=False
    )
    print("Step 3: Repression simulation complete.")
    print(delays_df)


if __name__ == "__main__":
    main()