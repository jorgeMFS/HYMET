#!/usr/bin/perl
use strict;
use warnings;
use Time::HiRes qw(time);

# User-modifiable parameters
my $mash_threshold_refseq = 0.90;  # Threshold inicial para RefSeq
my $mash_threshold_gtdb = 0.90;    # Threshold inicial para GTDB
my $mash_threshold_custom = 0.90;  # Threshold inicial para a nova base de dados

my $classification_processes = 8;  # Default value: 4
# Defining the base path
my $base_path = '/mnt/storagelv/home/inesbrancomartins/Tese/tool1';

# Defining paths for the scripts
my $mash_script = "$base_path/scripts/mash.sh";
my $download_script = "$base_path/scripts/downloadDB.py";
my $retry_download_script = "$base_path/scripts/retry_download.py"; # tenho de apagar esta linha 
my $minimap_script = "$base_path/scripts/minimap2.sh";
my $classification_script = "$base_path/scripts/classification.py";
my $cleandf_script = "$base_path/scripts/cleandf.py";

# Defining common paths
my $input_dir = "/mnt/storagelv/home/inesbrancomartins/Tese/stateoftheart/database/refdb/fungi/outputGCFfiltrado/";
my $output_dir = "$base_path/output";
my $data_dir = "$base_path/data";

# Specific paths for each script
my %paths = (
    mash => {
        mash_screen => "$data_dir/sketch1.msh",
        gtdb_mash_screen => "$data_dir/sketch2.msh",
        custom_mash_screen => "$data_dir/sketch3.msh",  # Novo caminho para a terceira base de dados
        screen_tab => "$output_dir/screen.tab",
        gtdb_screen_tab => "$output_dir/gtdb_screen.tab",
        custom_screen_tab => "$output_dir/custom_screen.tab",  # Novo caminho para o arquivo de saída da terceira base de dados
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
    minimap => {
        reference_set => "$data_dir/downloaded_genomes/combined_genomes.fasta",
        nt_mmi => "$output_dir/reference.mmi",
        resultados_paf => "$output_dir/resultados.paf",
    },
    classification => {
        paf_file => "$output_dir/resultados.paf",
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

# Mash with RefSeq dataset
run_command("$mash_script $input_dir $paths{mash}{mash_screen} $paths{mash}{screen_tab} $paths{mash}{filtered_screen} $paths{mash}{sorted_screen} $paths{mash}{top_hits} $paths{mash}{selected_genomes} $mash_threshold_refseq");

# Mash with GTDB dataset
run_command("$mash_script $input_dir $paths{mash}{gtdb_mash_screen} $paths{mash}{gtdb_screen_tab} /tmp/gtdb_filtered.tab /tmp/gtdb_sorted.tab /tmp/gtdb_top_hits.tab /tmp/gtdb_selected_genomes.txt $mash_threshold_gtdb");
run_command("cat /tmp/gtdb_selected_genomes.txt >> $paths{mash}{selected_genomes}");
run_command("sort -u -o $paths{mash}{selected_genomes} $paths{mash}{selected_genomes}");

# Mash with Custom dataset
run_command("$mash_script $input_dir $paths{mash}{custom_mash_screen} $paths{mash}{custom_screen_tab} /tmp/custom_filtered.tab /tmp/custom_sorted.tab /tmp/custom_top_hits.tab /tmp/custom_selected_genomes.txt $mash_threshold_custom");
run_command("cat /tmp/custom_selected_genomes.txt >> $paths{mash}{selected_genomes}");
run_command("sort -u -o $paths{mash}{selected_genomes} $paths{mash}{selected_genomes}");

# Download genomes
run_command("python3 $download_script $paths{download}{genomes_file} $paths{download}{downloaded_genomes} $paths{download}{taxonomy_file} $paths{download}{cache_dir}");

run_command("$minimap_script $input_dir $paths{minimap}{reference_set} $paths{minimap}{nt_mmi} $paths{minimap}{resultados_paf}");

# Classification
run_command("python3 $classification_script --paf $paths{classification}{paf_file} --taxonomy $paths{classification}{taxonomy_file} --hierarchy $paths{classification}{hierarchy_file} --output $paths{classification}{output_file} --processes $classification_processes");

my $end_time = time();
my $execution_time = $end_time - $start_time;

print "All scripts were executed successfully.\n";
print "\nExecution Summary:\n";
print "Total execution time: ".sprintf("%.2f", $execution_time)." seconds\n";

