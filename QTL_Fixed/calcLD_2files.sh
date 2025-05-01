



# Assign input files and MT value to variables
burnin_file=$1
Selection_file=$2
mt_value=$3

# Function to process a file through the pipeline
process_file() {
    local input_file=$1
    local output_prefix=$2
    local mt_value=$3

    # Step 1: Filter the VCF file using the provided MT value
    bcftools filter -e"MT = $mt_value" "$input_file" -o "${output_prefix}_filtered.vcf"

    # Step 2: Calculate LD using plink2
    plink2 --vcf "${output_prefix}_filtered.vcf" --r2-phased cols=+dprime,+freq --ld-window 999999 --ld-window-kb 999999 --ld-window-r2 0 --out "${output_prefix}_ld_results"
}

# Process the burnin file with the specified MT value
process_file "$burnin_file" "burnin" "$mt_value"

# Process the second file with the specified MT value
process_file "$Selection_file" "Selection_file" "$mt_value"

Rscript ~/vcftesting/LDheatmap_2files.R "${burnin_file}_ld_results.vcor" "${Selection_file}_ld_results.vcor"

echo "Processing complete."