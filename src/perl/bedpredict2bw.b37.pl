my $p=$ARGV[0] or die;
my $input=$ARGV[1] or die;
my $genome=$ARGV[2] || "/jet/home/dnaase/shared/data/genomes/hg19/hg19.chrom.sizes";
my $methy="$p.methy.bedgraph";
my $cov="$p.cov.bedgraph";
my $m_count="$p.methy_count.bedgraph";


	if($input=~/\.gz$/){
		`zcat $input | grep -v "start" | sed 's/^chrM/MT/g;s/^chr//g;'| cut -f1-3,5 | bedtools sort -i stdin > $m_count`;
		`zcat $input | grep -v "start" | sed 's/^chrM/MT/g;s/^chr//g;'| cut -f1-3,6 | bedtools sort -i stdin > $cov`;
		`zcat $input | grep -v "start" | sed 's/^chrM/MT/g;s/^chr//g;'| cut -f1-4 | bedtools sort -i stdin > $methy`;
	}else{
`grep -v "start" $input | sed 's/^chrM/MT/g;s/^chr//g;'| cut -f1-3,5 | bedtools sort -i stdin > $m_count`;
`grep -v "start" $input | sed 's/^chrM/MT/g;s/^chr//g;'| cut -f1-3,6 | bedtools sort -i stdin > $cov`;
`grep -v "start" $input | sed 's/^chrM/MT/g;s/^chr//g;'| cut -f1-4 | bedtools sort -i stdin > $methy`;
	}

my $methy_bw="$p.methy.b37.bw";
my $cov_bw="$p.cov.b37.bw";
my $m_count_bw="$p.methy_count.b37.bw";

`bedGraphToBigWig $methy $genome $methy_bw`;
`bedGraphToBigWig $cov $genome $cov_bw`;
`bedGraphToBigWig $m_count $genome $m_count_bw`;

`unlink $methy`;
`unlink $cov`;
`unlink $m_count`;



