#!/usr/bin/perl
use strict;
use warnings;
use Cwd qw(abs_path);
use FindBin qw($Bin);
use File::Path qw(make_path);
use File::Spec::Functions qw(catdir catfile rel2abs);

sub find_resource_base {
    my @candidates;

    push @candidates, $ENV{HYMET_RESOURCES} if $ENV{HYMET_RESOURCES};

    my $bin_path = abs_path($Bin);
    push @candidates, $bin_path if defined $bin_path;

    if (defined $bin_path) {
        my $prefix = abs_path(catdir($bin_path, '..'));
        if (defined $prefix) {
            push @candidates, catdir($prefix, 'share', 'hymet');
            push @candidates, catdir($prefix, 'share', 'HYMET');
        }
    }

    for my $candidate (@candidates) {
        next unless defined $candidate;
        my $scripts_path = catdir($candidate, 'scripts');
        my $taxonomy_script = catfile($scripts_path, 'taxonomy_hierarchy.py');
        if (-f $taxonomy_script) {
            return ($candidate, $scripts_path);
        }
    }

    die "Unable to locate HYMET resources. Set the HYMET_RESOURCES environment variable to the directory containing the scripts folder.\n";
}

my ($resource_base, $scripts_dir) = find_resource_base();

my $home = $ENV{HOME} // '.';
my $work_dir = $ENV{HYMET_WORKDIR} // catdir($home, '.hymet');
$work_dir = rel2abs($work_dir);

my $taxonomy_files_dir = catdir($work_dir, 'taxonomy_files');
my $data_dir = catdir($work_dir, 'data');
my $output_dir = catdir($work_dir, 'output');
my $cache_dir = catdir($work_dir, 'cache');

make_path($taxonomy_files_dir, $data_dir, $output_dir, $cache_dir);

print "Using work directory: $work_dir\n";

my $taxdmp_url = "ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdmp.zip";

print "Downloading taxonomy files...\n";
my $taxdmp_zip = catfile($taxonomy_files_dir, 'taxdmp.zip');
system('wget', '-q', '-O', $taxdmp_zip, $taxdmp_url);

print "Unzipping taxonomy files...\n";
system('unzip', '-q', $taxdmp_zip, '-d', $taxonomy_files_dir);

print "Executing Python script...\n";
my $taxonomy_script = catfile($scripts_dir, 'taxonomy_hierarchy.py');
system('python3', $taxonomy_script, $taxonomy_files_dir, $data_dir);

print "Configuration completed.\n";