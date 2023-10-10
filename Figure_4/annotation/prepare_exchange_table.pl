# this script prepares a lookupt table for ENSMUSG IDs to SYMBOL IDs
# result is printed to STDOUT

my @t_names;

my %was;

# input GTF is expected in the same folder as this script 
open FILE, "gencode.vM27.annotation.gtf";

while (<FILE>) {
	@parts = split(/\t/, $_);
	

	#transcript annotations contain the wanted information
	if ($parts[2] eq "gene") {
		#split at one or more symbols
		@info_parts = split(/["]/, $parts[8]);

		$ensmust = "";
		$ensmusg = "";
		$symbol = "";

		#extract ENSMUST, ENSMUG, SYMBOL and chromosome information
		for ($i=0; $i<=$#info_parts; $i++) {
			if ($info_parts[$i] =~ /ENSMUSG/) {$ensmusg = $info_parts[$i]}
			if ($info_parts[$i] =~ /gene_name/) {$symbol = $info_parts[$i+1]}
		}
		
		#prepare a list of ENSMUSG for output
		if (not defined($was{$ensmusg})) {

			print "$ensmusg\t$symbol\t$parts[0]\n";
			
			$was{$ensmusg}=1;
		}
	}
}
close(FILE);

