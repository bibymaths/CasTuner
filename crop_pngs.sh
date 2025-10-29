#!/usr/bin/env bash
set -euo pipefail

in_dir="R_results/plots"
crop_dir="R_results/plots_cropped"
png_dir="R_results/pngs"

mkdir -p "$crop_dir" "$png_dir"

# Check tools
command -v pdfcrop >/dev/null 2>&1 || { echo "ERROR: pdfcrop not found. Install texlive-extra-utils."; exit 1; }

if command -v magick >/dev/null 2>&1; then
  IM="magick"
elif command -v convert >/dev/null 2>&1; then
  IM="convert"
else
  echo "ERROR: ImageMagick not found. Install imagemagick."
  exit 1
fi

echo "Cropping PDFs in ${in_dir} ..."
shopt -s nullglob
for f in "$in_dir"/*.pdf; do
  base="$(basename "$f")"
  out_pdf="${crop_dir}/${base}"
  echo "  pdfcrop: $base"
  pdfcrop "$f" "$out_pdf"
done

echo "Converting cropped PDFs to PNG in ${png_dir} ..."
for f in "$crop_dir"/*.pdf; do
  base_png="${png_dir}/$(basename "${f%.pdf}.png")"
  echo "  rasterize: $(basename "$f") -> $(basename "$base_png")"
  $IM -density 200 "$f" -quality 100 png:- | $IM - -trim +repage "$base_png"
done

echo "Done. PNGs ready at: $png_dir"
