#!/bin/bash
## vizbin_count_table.sh


##### Script info

## vizbin_count_table: Generate count_table of exported vizbin sub-contigs
# Identifies which vizbin cluster the sub-contig is contained in
# Also counts how many times any sub-contigs from the parent contig occur in each vizbin cluster
# Use this table to select which bins to excluded from the filtering step that follows (i.e. which ones to keep in the data set)

## Required input 

# -i path/to/vizbin/out/files : path to output files from vizbin selection process (incl. contigs_[1-n].fna files and cluster_[1-n].fna files of whole clusters). (NOTE: the default file extension expected is .fna)

## Optional inputs

# -s vb_subcontig_export_prefix : prefix for files of selected problematic sub-contigs to putatively omit from bins (default = "contig")
# -c vb_cluster_prefix : prefix for files of whole vizbin clusters exported from vizbin (default = "cluster") 
# -e vb_files_extension: file extension for fasta files exported from VizBin for both cluster and contigs files (default = fna)


##### Set input options

# Set defaults for input fields
inpath=""
contig_prefix="contig"
cluster_prefix="cluster"
extension="fna"

# Help/usage documentation
function usage {
    echo -e \\n"Documentation for vizbin_count_table.sh"\\n
    echo -e "Required argument"\\n
    echo "-i path/to/vizbin/exported/files : path to output files from vizbin selection process (cluster and contigs files)."
    echo -e \\n"Optional arguments"\\n
    echo "-s contigs_prefix : prefix for files of selected problematic sub-contigs to putatively omit from bins (default = contig)."
    echo "-c clusters_prefix : prefix for files of whole vizbin clusters exported from vizbin (default = "cluster")."
    echo "-e file_extension : file extension for fasta files exported from VizBin for both cluster and contigs files (default = fna)."
    echo -e "-h : Display help message and exit."\\n
    echo -e "Example: ./vizbin_count_table.sh -i path/to/vizbin/exported/files/ -s contig -c cluster -e fna"\\n
    exit 1
}

# Incorporate commandline options
while getopts i:s:c:e:h option; do
    case "${option}" in
        i) inpath=${OPTARG%/};;
        s) contig_prefix=${OPTARG};;
        c) cluster_prefix=${OPTARG};;
        e) extension=${OPTARG};;
        h | [?]) usage ; exit;;
    esac
done


##### Error handling

# Check that inpath has been provided
if [ -z "${inpath}" ]; then
    echo "Error: path to input files not found"
    echo "< -i path/to/vizbin/out/files > must be provided to vizbin_count_table.sh"
    exit 1 
fi

# Check files matching cluster_prefix are present in the inpath directory
if ls ${inpath}/${contig_prefix}*.${extension} 1> /dev/null 2>&1; then
    :
else
    echo "Error: no input files found matching the provided file prefix (and/or file extension) for exported selected vizbin sub-contigs (defaults: prefix = selection; extension = fna)"
    exit 1
fi

# Check files matching cluster_prefix are present in the inpath directory
if ls ${inpath}/${cluster_prefix}*.${extension} 1> /dev/null 2>&1; then
    :
else
    echo "Error: no input files found matching the provided file prefix (and/or file extension) for exported vizbin clusters (defaults: prefix = cluster; extension = fna)"
    exit 1
fi


##### Main script    
    
# Create working list of parent contig IDs based on sub-contigs exported from vizbin
## 1. Create single fasta file of identified sub-contigs from multiple exports
## 2. extract sequence headers
## 3. trim sub-contig trailing numbers to retain only parent contig IDs
## 4. sort and retain only unique entries
cat ${inpath}/${contig_prefix}*.${extension} | grep ">" | sed -E -e "s/(cov_[0-9]+\.[0-9]+)\..*/\1/g" -e "s/>//g" | sort -u > vb_omit_contigs_tmp.txt

# Set up headers of vb_count_table.txt
echo -e -n "Subcontig_ID\tSubcontig_vb_cluster" > vb_count_table.txt
for vb_cluster in ${inpath}/${cluster_prefix}*.${extension}; do
    clustername=$(basename ${vb_cluster} .${extension})
    echo -e -n "\t${clustername}_count" >> vb_count_table.txt
done
echo -e "\tTotal_count" >> vb_count_table.txt

# loop over files containing extracted sub-contigs from vizbin (contigs_[1-n].fna files)
# For each:
## identify which vb_cluster it's located in (cluster_[1-n].fna)
## count how many times the parent contig of each sub-contig appears in each of the vb_clusters
## output this to vb_count_table.txt
while read -r line ; do
    echo "Processing ${line}"
    # Add sub-contig ID to count_table
    echo -e -n "${line}" >> vb_count_table.txt
    # Identify which vb_cluster file the sub-contig is in
    for vb_cluster in ${inpath}/${cluster_prefix}*.${extension}; do
        if grep -Fxq "${line}" ${vb_cluster}
        then
            # add cluster ID to count_table
            echo -e -n "\t$(basename ${vb_cluster} .${extension})" >> vb_count_table.txt
        else
            continue
        fi
    done
    # Strip sub-contig header back to parent contig ID
    parent=$(echo ${line} | sed -E "s/(cov_[0-9]+\.[0-9]+)\..*/\1/g")
    # count number of times it occurs in each vb_cluster
    total_count=0
    for vb_cluster in ${inpath}/${cluster_prefix}*.${extension}; do
        # count occurences in file
        cluster_count=$(grep -c ${parent} ${vb_cluster})
        # Add to running total_count for this contig
        total_count=$((${total_count}+${cluster_count}))
        # Add count to count_table
        echo -e -n "\t${cluster_count}" >> vb_count_table.txt
    done
    # Add total_count for this contig, and move to next line in the while loop
    echo -e "\t${total_count}" >> vb_count_table.txt 
done < <(cat ${inpath}/${contig_prefix}*.${extension} | grep ">")


##### END
