my $input=$ARGV[0];
my $p=$ARGV[1];
my $cpu=$ARGV[2] || 4;

`samtools sort -n -T $p --threads $cpu -n $input | samtools view --threads $cpu -f 3 -F 3852 -bh - | bedtools bamtobed -bedpe -mate1 -i stdin | perl -ne 'chomp;\@f=split "\t";if(\$f[0] ne \$f[3]){next;}\$s=\$f[1];\$e=\$f[5];if(\$f[8] eq "-"){\$s=\$f[4];\$e=\$f[2];}if(\$e>\$s){print "\$f[0]\\t\$s\\t\$e\\t.\\t\$f[7]\\t\$f[8]\\n";}' | sort -k1,1V -k2,2g | bgzip -fc --threads $cpu > $p.tsv.gz && tabix -f -p bed -b 2 -e 3 $p.tsv.gz`;
