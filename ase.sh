# [ase] alignment score extractor
find . -type f -path '*_output/*.txt' -exec grep -i "alignment score" {} + | awk -F": " '{print $2}'