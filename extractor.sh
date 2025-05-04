# [extractor] alignment score
find . -type f -path '*_output/*.txt' -exec grep -i "alignment score" {} + | awk -F": " '{print $2}'

# [extractor] E-value
find . -type f -path '*_output/*.txt' -exec grep -i "estimated e-value" {} + | awk -F": " '{print $2}'