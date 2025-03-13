#!/usr/bin/perl
use strict;
use warnings;
use Time::HiRes qw(time);

# User-modifiable parameters
my $mash_threshold_refseq = 0.90;  # Threshold for RefSeq
my $mash_threshold_gtdb = 0.83;   # Threshold for GTDB

my $classification_processes = 8;  # Default value: 4

# Defining the base path
my $base_path = '/mnt/storagelv/home/inesbrancomartins/Tese/tool1_fast';

# Defining paths for the scripts
my $mash_script = "$base_path/scripts/mash.sh";
my $download_script = "$base_path/scripts/downloadDB.py";
my $mashmap_script = "$base_path/scripts/mashmap.sh";  
my $classification_script = "$base_path/scripts/classification.py";
my $cleandf_script = "$base_path/scripts/cleandf.py";

# Defining common paths
my $input_dir = "/mnt/storagelv/home/inesbrancomartins/Tese/stateoftheart/database/refdb/vertebrateother/outputGCFfiltrado5";
my $output_dir = "$base_path/output";
my $data_dir = "$base_path/data";

# Specific paths for each script
my %paths = (
    mash => {
        mash_screen => "$data_dir/RefSeq88n.msh",
        gtdb_mash_screen => "$data_dir/sketcheddb.msh",
        screen_tab => "$output_dir/screen.tab",
        gtdb_screen_tab => "$output_dir/gtdb_screen.tab",
        filtered_screen => "$output_dir/filtered_screen.tab",
        sorted_screen => "$output_dir/sorted_screen.tab",
        top_hits => "$output_dir/top_hits.tab",
        selected_genomes => "$output_dir/selected_genomes.txt",
    },
    download => {
        genomes_file => "$output_dir/selected_genomes.txt",
        downloaded_genomes => "$data_dir/downloaded_genomes",
        taxonomy_file => "$data_dir/detailed_taxonomy.tsv",
        cache_dir => "$base_path/cache",
        log_file => "$data_dir/downloaded_genomes/genome_download_log.txt",
    },
    mashmap => {  
        reference_set => "$data_dir/downloaded_genomes/combined_genomes.fasta",
        resultados_mashmap => "$output_dir/mashmap.out",  # MashMap output
    },
    classification => {
        paf_file => "$output_dir/mashmap.out",  # Now uses MashMap output
        taxonomy_file => "$data_dir/detailed_taxonomy.tsv",
        hierarchy_file => "$data_dir/taxonomy_hierarchy.tsv",
        output_file => "$output_dir/classified_sequences.tsv",
    },
);

sub run_command {
    my ($command) = @_;
    print "Executing: $command\n";
    my $output = `$command 2>&1`;
    my $exit_status = $? >> 8;
   
    if ($exit_status == 0) {
        print "Success: $output\n";
    } else {
        print "Error executing: $output\n";
        exit($exit_status);
    }
}

my $start_time = time();

# Step 1: Run Mash with RefSeq dataset
print "Running Mash with RefSeq sketched database...\n";
run_command("$mash_script $input_dir $paths{mash}{mash_screen} $paths{mash}{screen_tab} $paths{mash}{filtered_screen} $paths{mash}{sorted_screen} $paths{mash}{top_hits} $paths{mash}{selected_genomes} $mash_threshold_refseq");

# Step 2: Run Mash with GTDB dataset
print "Running Mash with the costum sketched dataset...\n";
run_command("$mash_script $input_dir $paths{mash}{gtdb_mash_screen} $paths{mash}{gtdb_screen_tab} /tmp/gtdb_filtered.tab /tmp/gtdb_sorted.tab /tmp/gtdb_top_hits.tab /tmp/gtdb_selected_genomes.txt $mash_threshold_gtdb");

# Step 3: Combine selected genomes from all databases
print "Combining selected genomes from all databases...\n";
run_command("cat /tmp/gtdb_selected_genomes.txt >> $paths{mash}{selected_genomes}");
run_command("sort -u -o $paths{mash}{selected_genomes} $paths{mash}{selected_genomes}");

# Step 4: Download genomes
print "Downloading genomes...\n";
run_command("python3 $download_script $paths{download}{genomes_file} $paths{download}{downloaded_genomes} $paths{download}{taxonomy_file} $paths{download}{cache_dir}");

# # # Step 5: Concatenate input genomes into a single .fasta file for MashMap
my $concatenated_input = "$output_dir/concatenated_input.fasta";
run_command("cat $input_dir/*.fna > $concatenated_input");
print "Input genomes concatenated into $concatenated_input for MashMap.\n";

# Step 6: Run MashMap
print "Running MashMap...\n";
run_command("$mashmap_script $concatenated_input $paths{mashmap}{reference_set} $paths{mashmap}{resultados_mashmap}");

# Step 7: Classification
print "Running classification...\n";
run_command("python3 $classification_script --mashmap $paths{classification}{paf_file} --taxonomy $paths{classification}{taxonomy_file} --hierarchy $paths{classification}{hierarchy_file} --output $paths{classification}{output_file} --processes $classification_processes");

my $end_time = time();
my $execution_time = $end_time - $start_time;

print "All scripts were executed successfully.\n";
print "\nExecution Summary:\n";
print "Total execution time: ".sprintf("%.2f", $execution_time)." seconds\n";

