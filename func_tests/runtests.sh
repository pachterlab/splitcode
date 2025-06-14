#!/bin/bash

### Setup ###

splitcode="./src/splitcode --assign"
test_dir="./func_tests"

cmdexec() {
	cmd="$1"
	expect_fail="$2"
	printf "$cmd\n"
	if eval "$cmd 2> /dev/null 1> /dev/null"; then
		if [ "$expect_fail" = "1" ]; then
			printf "^[Failed - Returned normally but expected failure]\n"
			exit 1
		else
			printf "^[OK]\n"
		fi
	else
		if [ "$expect_fail" = "1" ]; then
			printf "^[OK - Exited non-zero as expected]\n"
		else
			printf "^[Failed]\n"
			exit 1
		fi
	fi
}

checkcmdoutput() {
	cmd="$1"
	correct_md5="$2"
	printf "$cmd\n"
	output_md5=$(eval "$cmd"|md5sum|awk '{ print $1 }')
	if [ "$output_md5" = "$correct_md5" ]; then
		printf "^[Output OK]\n"
	else
		printf "^[Output incorrect! Expected: "
		printf "$correct_md5"
		printf " Actual: "
		printf "$output_md5"
		printf "]\n"
		exit 1
	fi	
}

# Test that program can be run

cmdexec "$splitcode --version"

# Test some nesting workflows

checkcmdoutput "$splitcode --trim-only -c $test_dir/nest_config.txt --pipe --mod-names $test_dir/nest.test.1.R1.fq $test_dir/nest.test.1.R1.fa" 20484a6b5ba1f3c10931417310925ca9
checkcmdoutput "$splitcode --trim-only -c $test_dir/nest_config_2.txt --pipe --out-fasta --mod-names -N 2 $test_dir/nest.test.1.R1.fq $test_dir/nest.test.1.R2.fq $test_dir/nest.test.2.R1.fq $test_dir/nest.test.2.R2.fq" 9f1f430f7db125657c52653b298ee138
checkcmdoutput "$splitcode --trim-only -s $test_dir/summary.nest.txt -c $test_dir/nest_config_3.txt --pipe --empty=NNN --out-fasta --mod-names -N 2 $test_dir/nest.test.1.R1.fq $test_dir/nest.test.1.R2.fq $test_dir/nest.test.2.R1.fq $test_dir/nest.test.2.R2.fq" 6dabc80d0b3d813d07ab68fd313d2772
checkcmdoutput "cat $test_dir/summary.nest.txt|head -10|wc -l|tr -d '\n'|tr -d '\r'|tr -d ' '" d3d9446802a44259755d38e6d163e820
checkcmdoutput "$splitcode --trim-only -c $test_dir/nest_config_2.txt -o $test_dir/out.nest.1.fq,$test_dir/out.nest.2.fq,$test_dir/out.nest.3.fq --mod-names -N 2 $test_dir/nest.test.1.R1.fq $test_dir/nest.test.1.R2.fq $test_dir/nest.test.2.R1.fq $test_dir/nest.test.2.R2.fq" d41d8cd98f00b204e9800998ecf8427e
checkcmdoutput "cat $test_dir/out.nest.1.fq" 46e08a3f066c9ae5a741a31cf6d785fb
checkcmdoutput "cat $test_dir/out.nest.2.fq" 547e739c16eaa5495d1e566146ca2b65
checkcmdoutput "cat $test_dir/out.nest.3.fq" 5da0ce90bdd9f45c71948fae908e12b1
checkcmdoutput "$splitcode --trim-only -c $test_dir/nest_config_3.txt --x-only --empty=NNN --out-fasta --mod-names -N 2 $test_dir/nest.test.1.R1.fq $test_dir/nest.test.1.R2.fq $test_dir/nest.test.2.R1.fq $test_dir/nest.test.2.R2.fq" d41d8cd98f00b204e9800998ecf8427e
checkcmdoutput "cat $test_dir/xx.fasta" 6dabc80d0b3d813d07ab68fd313d2772
checkcmdoutput "$splitcode --outb=$test_dir/bc -u $test_dir/a,$test_dir/b,$test_dir/c -o $test_dir/aa,$test_dir/bb,$test_dir/cc -m $test_dir/nest.mapping.txt -c $test_dir/nest_config_2.txt --out-fasta --mod-names --com-names -N 2 $test_dir/nest.test.1.R1.fq $test_dir/nest.test.1.R2.fq $test_dir/nest.test.2.R1.fq $test_dir/nest.test.2.R2.fq" d41d8cd98f00b204e9800998ecf8427e
checkcmdoutput "cat $test_dir/a" c1c334bec88086eae1cafe58857e286b
checkcmdoutput "cat $test_dir/b" 19f5e79683fdb9e96e3cc0981e13c0cb
checkcmdoutput "cat $test_dir/c" 66d3a7966397c33fa9bc0c3e7f69ebda
checkcmdoutput "cat $test_dir/aa" 0d81ab04a3ffef6cb904311f08f61e32
checkcmdoutput "cat $test_dir/bb" edf67aa6dab16e0953fbd9895285fb8b
checkcmdoutput "cat $test_dir/cc" 732e911683d3e5cd8ff2aad7fd1fafee
checkcmdoutput "cat $test_dir/bc" e39b5fc29e0dc621ac15c24262fcf507
checkcmdoutput "cat $test_dir/nest.mapping.txt" 562376738ccdb0645824021a1d91954b
checkcmdoutput "$splitcode --no-outb --pipe -u $test_dir/a,$test_dir/b,$test_dir/c -m $test_dir/nest.mapping.txt -c $test_dir/nest_config_2.txt --out-fasta --mod-names -N 2 $test_dir/nest.test.1.R1.fq $test_dir/nest.test.1.R2.fq $test_dir/nest.test.2.R1.fq $test_dir/nest.test.2.R2.fq" 0634a9a8aeff98e586e0dbb91a0e8f04
checkcmdoutput "cat $test_dir/b" 19f5e79683fdb9e96e3cc0981e13c0cb
checkcmdoutput "cat $test_dir/nest.mapping.txt" 562376738ccdb0645824021a1d91954b
checkcmdoutput "$splitcode -m $test_dir/nest.mapping.txt --nFastqs=2 -c $test_dir/nest_config_4.txt --mod-names --empty=N --outb=$test_dir/bc --x-only -u $test_dir/a,$test_dir/b $test_dir/nest.test.3.R1.fq $test_dir/nest.test.3.R2.fq" d41d8cd98f00b204e9800998ecf8427e
checkcmdoutput "cat $test_dir/bc" b4ef4481277e4aff5adaf3513c988e2b
checkcmdoutput "cat $test_dir/a" f501b402785d361b968b7ef66aa6100e
checkcmdoutput "cat $test_dir/b" 0f0a86618bd4b2fcf3996114ab437d42
checkcmdoutput "$splitcode --trim-only --no-x-out -c $test_dir/nest_config_5.txt --pipe --empty=NNN --out-fasta --mod-names -N 2 $test_dir/nest.test.1.R1.fq $test_dir/nest.test.1.R2.fq $test_dir/nest.test.2.R1.fq $test_dir/nest.test.2.R2.fq" b5138c11c6250d138252f81925ad3778
checkcmdoutput "$splitcode --trim-only -c $test_dir/nest_config_5.txt  --pipe --empty=NNN --out-fasta --mod-names -N 2 $test_dir/nest.test.1.R1.fq $test_dir/nest.test.1.R2.fq $test_dir/nest.test.2.R1.fq $test_dir/nest.test.2.R2.fq" 00667a64eb00f8b6993e2f2819e138d7



# Test some --barcode-encode

checkcmdoutput "$splitcode -B "B,D,E" -c $test_dir/config_fastassign.txt --pipe --mod-names -u $test_dir/unassign.fq -m $test_dir/mx.txt $test_dir/nest.test.3.R1.fq" 4e9181c048e722c8fb69fb8ac58f726e
checkcmdoutput "cat $test_dir/unassign.fq" 0eb7dccdeac130ea6c33e5174a7d189a
checkcmdoutput "cat $test_dir/mx.txt" 129a22a29e2ee9ae3b29389a20ed0bff

checkcmdoutput "$splitcode -B "B,D,E" --bclen=20 -c $test_dir/config_fastassign_2.txt --pipe --mod-names -u $test_dir/unassign.fq -m $test_dir/mx.txt $test_dir/nest.test.3.R1.fq" d3c5c22523d87e22ffdf11bf58adccba
checkcmdoutput "cat $test_dir/unassign.fq" d41d8cd98f00b204e9800998ecf8427e
checkcmdoutput "cat $test_dir/mx.txt" b95e8b332c8a0a7ffc0f91118e754302



# Test lift workflow

checkcmdoutput "$splitcode --lift $test_dir/vcf_validation.fa.gz $test_dir/test_1.vcf.gz CAST_EiJ --kmer-length=31 --kmer-output=$test_dir/test.kmers.1.txt" 7e9c1d67efdf7113bfab367cb5d2d640
checkcmdoutput "$splitcode --lift $test_dir/vcf_validation.fa.gz $test_dir/test_2.vcf.gz CAST_EiJ --rename --kmer-length=31 --kmer-output=$test_dir/test.kmers.2.txt" c9a91833da19b20383d6bd1d3e32ff8f
checkcmdoutput "cat $test_dir/test.kmers.1.txt" 54d189f4549f6b35ea80ec5c167332b7
checkcmdoutput "cat $test_dir/test.kmers.2.txt" 9982de087ed358724580836600cc9ba7
checkcmdoutput "$splitcode --lift $test_dir/vcf_validation.fa.gz $test_dir/test_2.vcf.gz CAST_EiJ --diploid --kmer-length=31 --kmer-output=$test_dir/test.kmers.2.txt" 6c33bd3ba7cd9aaefeca5cbaa272cdfe
checkcmdoutput "cat $test_dir/test.kmers.2.txt" a4c73b67a0ad6e5094ff7f2dfdda15bc
checkcmdoutput "$splitcode --lift --snv-only $test_dir/vcf_validation.fa.gz $test_dir/test_2.vcf.gz CAST_EiJ" 4bdcaf9f34da45033d477651d8845bf2
checkcmdoutput "$splitcode --lift --kmer-sj $test_dir/vcf_validation.fa.gz $test_dir/example.SJ.tab --kmer-length=31 --kmer-header=X_ --kmer-header-num" 1e88d0b72323a6c11faa73e19007fc5f


# Test from-name, random, and revcomp

checkcmdoutput "$splitcode --trim-only --pipe --from-name=0,0,2,:: $test_dir/from_name.fq" e0b141cd0e6348d3cd64d23d61f67d09
checkcmdoutput "$splitcode --trim-only --pipe --from-name=\"0,0,0,::;0,0,0,::+\" -b GTTTA -d 1 --revcomp=1 --mod-names $test_dir/from_name.fq" 83b2f368aea4e889ad5995dbaef6f744
checkcmdoutput "$splitcode --trim-only --pipe --random=0,2,25 $test_dir/from_name.fq" 669d8dbda3e81f74159e1ce280686d49
checkcmdoutput "$splitcode --trim-only --pipe  -x \"<umi{:}>\" --random=0,-1,25 $test_dir/from_name.fq" 30b84f54f195c1e499cafb33acb1b2bf
checkcmdoutput "$splitcode --trim-only --pipe --from-name=\"0,0,-1,::;0,0,-1,::+\" -x \"<umi{:}>,<umi{:}>\" $test_dir/from_name.fq" 21f22aeb47e01c059f957fddb28a01aa
checkcmdoutput "$splitcode --trim-only --pipe --from-name=\"0,0,-1,::;0,0,-1,::+\" -x \"<umi{:}>,<umi2{:}>\" $test_dir/from_name.fq" c5456433bca2f5599d15f0785637149e



# Test SPRITE config files

cmdexec "$splitcode -t 1 -o "$test_dir/A_out_1.fastq.gz,$test_dir/A_out_2.fastq.gz" -O $test_dir/A_out_barcodes.fastq.gz -N 2 -c $test_dir/splitcode_example_config.txt -y <(echo "DPM,Y,ODD,EVEN,ODD") --mod-names --gzip -m $test_dir/A_out_mapping.txt.gz $test_dir/A_1.fastq.gz $test_dir/A_2.fastq.gz"
checkcmdoutput "zcat < $test_dir/A_out_1.fastq.gz" 05eb0a03078777fb36ec0a463707b1c6
checkcmdoutput "zcat < $test_dir/A_out_2.fastq.gz" 6780951145cb17370b8977aad4513355
checkcmdoutput "zcat < $test_dir/A_out_barcodes.fastq.gz" 395253488fa5a04acdb42af9291e1139
checkcmdoutput "zcat < $test_dir/A_out_mapping.txt.gz" 3bed011f0405480768019a0f2ceb445e
checkcmdoutput "$splitcode -t 1 --pipe -N 2 -c $test_dir/splitcode_example_config.txt -y <(echo "DPM,Y,ODD,EVEN,ODD") --mod-names -m /dev/null $test_dir/A_1.fastq.gz $test_dir/A_2.fastq.gz" a9857b5416fa815d63f3b262252f5f4f
checkcmdoutput "$splitcode -t 1 --pipe -N 2 -c $test_dir/splitcode_example_config.txt -y <(echo "DPM,Y,ODD,EVEN,ODD") --com-names -m /dev/null $test_dir/A_1.fastq.gz $test_dir/A_2.fastq.gz" d9bd0b3183cfb50256d5473af2f6476e
checkcmdoutput "$splitcode -t 1 -y <(echo "Y,ODD,EVEN,LIGTAG") -N 2 -c $test_dir/splitcode_example_config_2.txt --mod-names --com-names --pipe -m /dev/null $test_dir/B_1.fastq.gz $test_dir/B_2.fastq.gz" af7415d5cae5cd87ec9c38faeb0a3093
checkcmdoutput "$splitcode -t 1 -y <(echo "Y,ODD,EVEN,LIGTAG") -N 2 -c $test_dir/splitcode_example_config_3.txt --mod-names --com-names --pipe -m /dev/null $test_dir/B_1.fastq.gz $test_dir/B_2.fastq.gz" af7415d5cae5cd87ec9c38faeb0a3093
checkcmdoutput "$splitcode -t 1 -y <(echo "Y,ODD,EVEN,LIGTAG") -N 2 -c $test_dir/splitcode_example_config_3.txt -P TTCC -n 10 --mod-names --pipe -m /dev/null $test_dir/B_1.fastq.gz $test_dir/B_2.fastq.gz" 964a52f977964e1c2eeb4ffadf49bdcd

# Testing files supplied in config file

checkcmdoutput "$splitcode -t 1 --no-outb -m /dev/null -c $test_dir/configb.txt --pipe --mod-names $test_dir/b.fastq" 86479a52fd7d35a1d0fcb62c5dbde424
cat $test_dir/splitcode_example_config_4.txt > $test_dir/splitcode_example_config_4_.txt
echo "EVEN"$'\t'"Even2Bo"$'\t'"${test_dir}/even.txt\$"$'\t'"2"$'\t'"1:30"$'\t'"1"$'\t'"{{ODD}}6-12
ODD"$'\t'"Odd2Bo"$'\t'"${test_dir}/odd.txt\$"$'\t'"2"$'\t'"1:15"$'\t'"2"$'\t'"{{EVEN}}6-12" >> $test_dir/splitcode_example_config_4_.txt
checkcmdoutput "$splitcode -t 1 --pipe -N 2 -c $test_dir/splitcode_example_config_4_.txt -y <(echo "DPM,Y,ODD,EVEN,ODD") -m /dev/null $test_dir/A_1.fastq.gz $test_dir/A_2.fastq.gz" 5ff23e700eeb1a86eab3c0cf60a39108
checkcmdoutput "$splitcode -t 1 --pipe -N 2 -c $test_dir/splitcode_example_config_4_.txt -y <(echo "DPM,Y,ODD,EVEN,ODD") -n 100 -m /dev/null $test_dir/A_1.fastq.gz $test_dir/A_2.fastq.gz" b2da5c0f4a2e3ade60583c31eb5ac93d
checkcmdoutput "$splitcode -t 1 --pipe -N 2 -c $test_dir/splitcode_example_config_4_.txt -y <(echo "DPM,Y,ODD,EVEN,ODD") -n 100 --mod-names -m /dev/null $test_dir/A_1.fastq.gz $test_dir/A_2.fastq.gz" 13b5021426756b65a39aea31f1518186
cat $test_dir/splitcode_example_config_4.txt > $test_dir/splitcode_example_config_4__.txt
echo "EVEN"$'\t'"Even2Bo"$'\t'"${test_dir}/even.txt"$'\t'"2"$'\t'"1:30"$'\t'"1"$'\t'"{{ODD}}6-12
ODD"$'\t'"Odd2Bo"$'\t'"${test_dir}/odd.txt"$'\t'"2"$'\t'"1:15"$'\t'"2"$'\t'"{{EVEN}}6-12" >> $test_dir/splitcode_example_config_4__.txt
checkcmdoutput "$splitcode -t 1 --pipe -N 2 -c $test_dir/splitcode_example_config_4__.txt -y <(echo "DPM,Y,ODD,EVEN,ODD") -n 100 -m /dev/null $test_dir/A_1.fastq.gz $test_dir/A_2.fastq.gz" b2da5c0f4a2e3ade60583c31eb5ac93d
checkcmdoutput "$splitcode -t 1 --pipe -N 2 -c $test_dir/splitcode_example_config_4__.txt -y <(echo "DPM,Y,ODD,EVEN,ODD") -n 100 --mod-names -m /dev/null $test_dir/A_1.fastq.gz $test_dir/A_2.fastq.gz" 3fad75cf39a30f328434370c70f995a0


# Basic UMI extraction testing

checkcmdoutput "$splitcode -t 1 -N 2 -c $test_dir/splitcode_example_config.txt -y <(echo "DPM,Y,ODD,EVEN,ODD") -m /dev/null --pipe -x \"{Odd2Bo50}<umi[6]>\" $test_dir/A_1.fastq.gz $test_dir/A_2.fastq.gz|head -16" 4a413abb788efb880a35e6d60ca4ac4d
checkcmdoutput "$splitcode -t 1 -N 2 -c $test_dir/splitcode_example_config.txt -y <(echo "DPM,Y,ODD,EVEN,ODD") -m /dev/null --pipe -x \"{Odd2Bo50}<~umi[6]>\" $test_dir/A_1.fastq.gz $test_dir/A_2.fastq.gz|head -16" a2c4930122a996c717d715080c8f3b22
checkcmdoutput "$splitcode -t 1 -N 2 -c $test_dir/splitcode_example_config.txt -y <(echo "DPM,Y,ODD,EVEN,ODD") --x-names -m /dev/null --pipe -x \"{Odd2Bo50}<umi[6]>\" $test_dir/A_1.fastq.gz $test_dir/A_2.fastq.gz|head -16" bea2940a59a5fcc833310e894a9b2c2d
checkcmdoutput "$splitcode -t 1 -N 2 -c $test_dir/splitcode_example_config.txt -y <(echo "DPM,Y,ODD,EVEN,ODD") --x-names -m /dev/null --pipe -x \"{Odd2Bo50}<umi[6]>\" --out-fasta --empty-remove $test_dir/A_1.fastq.gz $test_dir/A_2.fastq.gz|head -16" 645eb8593b524db799d096eaa60c165e

# Advanced UMI extraction testing 1

checkcmdoutput "$splitcode -t 1 -N 2 -c $test_dir/splitcode_example_config.txt -y <(echo "DPM,Y,ODD,EVEN,ODD") --x-names -m /dev/null --pipe -x \"{Odd2Bo50}<umi[6-120]>,{Odd2Bo50}<umi[10-100]>,{Odd2Bo50}<umi[6-9]>\" $test_dir/A_1.fastq.gz $test_dir/A_2.fastq.gz|head -16" 9c78c856d1bb693b784cf44c6acea43f
checkcmdoutput "$splitcode -t 1 -N 2 -c $test_dir/splitcode_example_config.txt -y <(echo "DPM,Y,ODD,EVEN,ODD") --x-names -m /dev/null --pipe -x \"{Odd2Bo50}<umi1[39]>,{Odd2Bo50}<umi2[40]>\" $test_dir/A_1.fastq.gz $test_dir/A_2.fastq.gz|head -20" 27ac575249c7338c3bf90ab1ab1f32c0
checkcmdoutput "$splitcode -t 1 -N 2 -c $test_dir/splitcode_example_config.txt -y <(echo "DPM,Y,ODD,EVEN,ODD") --x-names -m /dev/null --pipe -x \"{Odd2Bo50}<umi[41]>\" $test_dir/A_1.fastq.gz $test_dir/A_2.fastq.gz|head -16" 03ba76e5c5f2649a6e39d71490b81f73
checkcmdoutput "$splitcode -t 1 -N 2 -c $test_dir/splitcode_example_config.txt -y <(echo "DPM,Y,ODD,EVEN,ODD") --x-names -m /dev/null --pipe -x \"{Odd2Bo50}<umi[6]>,{Odd2Bo50}<umi[9]>,{Odd2Bo50}9<umi1[5]>,{Odd2Bo50}<umi2[9]>,{Odd2Bo50}<umi1[3]>\" $test_dir/A_1.fastq.gz $test_dir/A_2.fastq.gz|head -24" 5481a84873d0fc929a28ba47ed6a9f6e
checkcmdoutput "$splitcode -t 1 -N 2 -c $test_dir/splitcode_example_config.txt -y <(echo "DPM,Y,ODD,EVEN,ODD") --x-names -m /dev/null --pipe -x \"{Odd2Bo50}<umi[6]>,{Odd2Bo50}<umi[9]>,{Odd2Bo50}9<umi1[5]>,{Odd2Bo50}<umi2[9]>,{Odd2Bo50}<umi1[3]>\" --no-chain $test_dir/A_1.fastq.gz $test_dir/A_2.fastq.gz|head -24" 76c026646e7de5b06ff14dbc4b142173
checkcmdoutput "$splitcode -t 1 -N 2 -c $test_dir/splitcode_example_config.txt -y <(echo "DPM,Y,ODD,EVEN,ODD") --no-outb --x-only -m /dev/null --pipe -x \"{Odd2Bo50}<umi[6]>,{Odd2Bo50}<umi[9]>,{Odd2Bo50}9<umi1[5]>,{Odd2Bo50}<umi2[9]>,{Odd2Bo50}<umi1[3]>\" $test_dir/A_1.fastq.gz $test_dir/A_2.fastq.gz|head -24" dbbad38ef87df8448481128232f879e6
checkcmdoutput "$splitcode -t 1 -N 2 -c $test_dir/splitcode_example_config.txt -y <(echo "DPM,Y,ODD,EVEN,ODD") --x-names -m /dev/null --pipe -x \"{Even2Bo33}<umi>{Odd2Bo50},{Even2Bo33}1<umi>2{Odd2Bo50},{Even2Bo33}1<umi[4-5]>2{Odd2Bo50},{Even2Bo33}1<umi[7-9]>2{Odd2Bo50},{Even2Bo33}1<umi1[5-7]>2{Odd2Bo50},{Even2Bo33}1<umi1[6]>2{Odd2Bo50},{Even2Bo33}1<umi2[6]>2{Odd2Bo50}\" $test_dir/A_1.fastq.gz $test_dir/A_2.fastq.gz|head -24" ef843b20cc6af4f692840692ffffe680
checkcmdoutput "$splitcode -t 1 -N 2 -c $test_dir/splitcode_example_config.txt -y <(echo "DPM,Y,ODD,EVEN,ODD") --x-names -m /dev/null --pipe -x \"{DPM6H4}1<umi0[5]>,{DPM6H4}1<umi1>{Even2Bo33},{Odd2Bo10}<umi2>{Even2Bo33},{Odd2Bo10}<umi3>{Odd2Bo50},{Odd2Bo10}2<umi4>{Even2Bo33},{Odd2Bo10}<umi4[1-2]>2{Odd2Bo50},{Odd2Bo50}<umi5[5-100]>,{Odd2Bo10}<umi4>2{Odd2Bo50}\" $test_dir/A_1.fastq.gz $test_dir/A_2.fastq.gz|head -36" 85e0e59877d3402f14f57dba424f96ed

# Advanced UMI extraction testing 2

checkcmdoutput "$splitcode -t 1 -N 2 -c $test_dir/splitcode_example_config.txt -y <(echo "DPM,Y,ODD,EVEN,ODD") --x-names -m /dev/null --pipe -x \"{{ODD}}<umi>{{ODD}},{{ODD}}<umi1>{Odd2Bo50},{{ODD}}<umi2>{{EVEN}}\" $test_dir/A_1.fastq.gz $test_dir/A_2.fastq.gz|head -24" 4e1e5cd270efdacbe9de77b568fb03da
checkcmdoutput "$splitcode -t 1 -N 2 -c $test_dir/splitcode_example_config.txt -y <(echo "DPM,Y,ODD,EVEN,ODD") --x-names -m /dev/null --pipe -x \"{{DPM}}<umi>{{EVEN}},{{EVEN}}<umi>{{EVEN}},{{EVEN}}<umi>{Odd2Bo10}\" $test_dir/A_1.fastq.gz $test_dir/A_2.fastq.gz|head -16" 03ba76e5c5f2649a6e39d71490b81f73
checkcmdoutput "$splitcode -t 1 -N 2 -c $test_dir/splitcode_example_config.txt -y <(echo "DPM,Y,ODD,EVEN,ODD") --x-names -m /dev/null --pipe -x \"{{EVEN}}<umi1>{Odd2Bo50},{{ODD}}<umi2>{{EVEN}},{{ODD}}<umi3>{{ODD}},{{EVEN}}<umi4>{{EVEN}},{{EVEN}}<umi5>{{ODD}}\" $test_dir/A_1.fastq.gz $test_dir/A_2.fastq.gz|head -32" c5f22e5f8d29c145a44bd95b0bbcd190

# Advanced UMI extraction testing 3

checkcmdoutput "$splitcode -t 1 -N 2 -c $test_dir/splitcode_example_config.txt -y <(echo "DPM,Y,ODD,EVEN,ODD") --x-names -m /dev/null --pipe -x \"0:4<umi[92]>,0:1<umi[10]>,0:1<umi1[5]>,0:3<umi2[5]>,0:4<umi3[6]>,1:1<umi4[5]>,0:5<umi5[8]>\" $test_dir/A_1.fastq.gz $test_dir/A_2.fastq.gz|head -36" b7f129379a232040c83a638961619644
checkcmdoutput "$splitcode -t 1 -N 2 -c $test_dir/splitcode_example_config.txt -y <(echo "DPM,Y,ODD,EVEN,ODD") --x-names -m /dev/null --pipe -x \"<umi[10]>0:1,<umi1[10]>0:10,<umi1[10]>0:20,<umi2[20]>0:-1,<umi3[20]>1:-1\" $test_dir/A_1.fastq.gz $test_dir/A_2.fastq.gz|head -28" 39aae648d41c03eaa3c2ccc4948866aa
checkcmdoutput "$splitcode -t 1 -N 2 -c $test_dir/splitcode_example_config.txt -y <(echo "DPM,Y,ODD,EVEN,ODD") --x-names -m /dev/null --pipe -x \"<umi1[96]>0:-1,<umi2[10]>0:10,<umi2[10]>0:10,<umi3[20]>1:-1,<umi4[21]>1:-1\" $test_dir/A_1.fastq.gz $test_dir/A_2.fastq.gz|head -28" 466ac29defbd770b78c4a0fdf25bd5a4
checkcmdoutput "$splitcode -t 1 -N 2 -c $test_dir/splitcode_example_config.txt -y <(echo "DPM,Y,ODD,EVEN,ODD") --x-names -m /dev/null --pipe -x \"{{EVEN}}<umi1>1:80,{{EVEN}}<umi2>1:-1\" $test_dir/A_1.fastq.gz $test_dir/A_2.fastq.gz|head -20" e42fe13d5265227af249325714495e75
checkcmdoutput "$splitcode -t 1 -N 2 -c $test_dir/splitcode_example_config.txt -y <(echo "DPM,Y,ODD,EVEN,ODD") --x-names -m /dev/null --pipe -x \"{{EVEN}}<umi1>1:1000\" $test_dir/A_1.fastq.gz $test_dir/A_2.fastq.gz|head -16" 03ba76e5c5f2649a6e39d71490b81f73
checkcmdoutput "$splitcode -t 1 -N 2 -c $test_dir/splitcode_example_config.txt -y <(echo "DPM,Y,ODD,EVEN,ODD") --x-names -m /dev/null --pipe -x \"{{ODD}}<umi1>1:-1\" $test_dir/A_1.fastq.gz $test_dir/A_2.fastq.gz|head -16" dcb672a30a80bc2b12dbf8b4a1df81f6
checkcmdoutput "$splitcode -t 1 -N 2 -c $test_dir/splitcode_example_config.txt -y <(echo "DPM,Y,ODD,EVEN,ODD") --x-names -m /dev/null --pipe -x \"{{EVEN}}5<umi1>1:-1\" $test_dir/A_1.fastq.gz $test_dir/A_2.fastq.gz|head -16" 2687f12d1eb72da8c2999ddbf45da38f
checkcmdoutput "$splitcode -t 1 -N 2 -c $test_dir/splitcode_example_config.txt -y <(echo "DPM,Y,ODD,EVEN,ODD") --x-names -m /dev/null --pipe -x \"{{EVEN}}5<umi1[60-70]>1:120,{{EVEN}}5<umi1[55-70]>1:120\" $test_dir/A_1.fastq.gz $test_dir/A_2.fastq.gz|head -16" 2687f12d1eb72da8c2999ddbf45da38f

# Advanced UMI extraction testing 4

checkcmdoutput "$splitcode -t 1 -N 2 -c $test_dir/splitcode_example_config.txt -y <(echo "DPM,Y,ODD,EVEN,ODD") --x-names -m /dev/null --pipe -x \"1:41<umi1>{{ODD}},1:18<umi2>{{ODD}},1:17<umi3>{{ODD}},1:16<umi4>{{ODD}}\" $test_dir/A_1.fastq.gz $test_dir/A_2.fastq.gz|head -28" e612346ccb3010380dfc93c9b02975ca
checkcmdoutput "$splitcode -t 1 -N 2 -c $test_dir/splitcode_example_config.txt -y <(echo "DPM,Y,ODD,EVEN,ODD") --x-names -m /dev/null --pipe -x \"1:15<umi1>{Odd2Bo10},1:15<umi2>{Odd2Bo50}\" $test_dir/A_1.fastq.gz $test_dir/A_2.fastq.gz|head -20" 07ae6220ccbd7099533d2b787d93ad1a
checkcmdoutput "$splitcode -t 1 -N 2 -c $test_dir/splitcode_example_config.txt -y <(echo "DPM,Y,ODD,EVEN,ODD") --x-names -m /dev/null --pipe -x \"1:15<umi1>1{Odd2Bo10},1:15<umi2>4{Odd2Bo50},1:15<umi3>4{Odd2Bo10}\" $test_dir/A_1.fastq.gz $test_dir/A_2.fastq.gz|head -24" 54e8f3ed651c4eeea9265010bfc58e31
checkcmdoutput "$splitcode -t 1 -N 2 -c $test_dir/splitcode_example_config.txt -y <(echo "DPM,Y,ODD,EVEN,ODD") --x-names -m /dev/null --pipe -x \"1:30<umi1>1{{ODD}},1:30<umi2>1{Odd2Bo10},1:1000<umi3>4{Odd2Bo50},0:30<umi4>{{ODD}}\" $test_dir/A_1.fastq.gz $test_dir/A_2.fastq.gz|head -28" 941b4e3bca1d84f61e19be27065c7094

# Advanced UMI extraction testing 5

checkcmdoutput "$splitcode -t 1 -N 2 -c $test_dir/splitcode_example_config.txt -y <(echo "DPM,Y,ODD,EVEN,ODD") --x-names -m /dev/null --pipe -x \"1:10<umi1>1:20,1:10<umi2[8-12]>1:20\" $test_dir/A_1.fastq.gz $test_dir/A_2.fastq.gz|head -20" 534996586a95f64431f6e9987e28307c
checkcmdoutput "$splitcode -t 1 -N 2 -c $test_dir/splitcode_example_config.txt -y <(echo "DPM,Y,ODD,EVEN,ODD") --x-names -m /dev/null --pipe -x \"1:10<umi1[5]>1:20,1:10<umi2[3-8]>1:20,1:10<umi3[12-16]>1:20\" $test_dir/A_1.fastq.gz $test_dir/A_2.fastq.gz|head -24" 6c5146a762f06f264c4a7e43bc87cc2e
checkcmdoutput "$splitcode -t 1 -N 2 -c $test_dir/splitcode_example_config.txt -y <(echo "DPM,Y,ODD,EVEN,ODD") --x-names -m /dev/null --pipe -x \"1:10<umi1[5]>1:20,1:10<umi2[3-8]>1:20,1:10<umi3[12-16]>1:20\" --empty-remove $test_dir/A_1.fastq.gz $test_dir/A_2.fastq.gz|head -24" fb989ecb89d6ec059d824116c12509bd
checkcmdoutput "$splitcode -t 1 -N 2 -c $test_dir/splitcode_example_config.txt -y <(echo "DPM,Y,ODD,EVEN,ODD") --x-names -m /dev/null --pipe -x \"1:10<umi1[8-12]>1:20,1:15<umi2[8-12]>1:25,1:15<umi3[12-20]>1:25,1:15<umi4[8-12]>1:25,1:15<umi4[8-12]>1:25\" $test_dir/A_1.fastq.gz $test_dir/A_2.fastq.gz|head -28" 3e6d524540806d03bcaa09ea3991406e

# Quality trimming tests

echo "@read0
AAGCTTCCGG
+
KI;<)(,%#$
@read1
AAGCTTCCGG
+
I;<*,(,%#$" > $test_dir/test.fq


checkcmdoutput "$splitcode --trim-only --pipe $test_dir/test.fq" cb52b79ed7469ca2ffe5739ec544b157
checkcmdoutput "$splitcode --trim-only --pipe -q 10 --qtrim-3 $test_dir/test.fq" 0284bb3b6fa601d6f85d8d051f7e6431
checkcmdoutput "cat $test_dir/test.fq|$splitcode --trim-only --pipe -q 10 --qtrim-3 -" 0284bb3b6fa601d6f85d8d051f7e6431
checkcmdoutput "$splitcode --trim-only --pipe -q 10 --qtrim-5 --qtrim-3 $test_dir/test.fq" 0284bb3b6fa601d6f85d8d051f7e6431
checkcmdoutput "$splitcode --trim-only --pipe -q 10 --qtrim-5 $test_dir/test.fq" cb52b79ed7469ca2ffe5739ec544b157
checkcmdoutput "$splitcode --trim-only --pipe -q 10 --qtrim-3 --qtrim-naive $test_dir/test.fq" 3b25d6c8a52ad921c4f6565ac7b1b299
checkcmdoutput "$splitcode --trim-only --pipe -q 11 --qtrim-3 --qtrim-naive $test_dir/test.fq" 3b25d6c8a52ad921c4f6565ac7b1b299
checkcmdoutput "$splitcode --trim-only --pipe -q 12 --qtrim-3 --qtrim-naive $test_dir/test.fq" 903366847aa3aa415250f57d83f5517d
checkcmdoutput "$splitcode --trim-only --pipe -q 10 --qtrim-5 --qtrim-naive $test_dir/test.fq" cb52b79ed7469ca2ffe5739ec544b157
checkcmdoutput "$splitcode --trim-only --pipe -q 8 --qtrim-3 --qtrim-naive --phred64 $test_dir/test.fq" f03fb0cee9154be006dfccd5f67797f6

# Quality trimming tests - Advanced

checkcmdoutput "$splitcode --trim-only --pipe -q 10 --qtrim-3 --qtrim-naive -x \"0:2<umi1>0:-1\" $test_dir/test.fq" 507cd03e7a97e95afc1ea1018976a623
checkcmdoutput "$splitcode --trim-only --pipe -q 10 --qtrim-3 --qtrim-naive --qtrim-pre -x \"0:2<umi1>0:-1\" $test_dir/test.fq" 021a446d088d3dfe37028873f223bc0a
checkcmdoutput "$splitcode --trim-only --pipe -q 10 --qtrim-3 --qtrim-pre -x \"0:2<umi1>0:-1\" $test_dir/test.fq" 9f8a7ebcbb903f5828f60a8c5d142bf3
checkcmdoutput "$splitcode --trim-only --pipe -q 10 --qtrim-3 --qtrim-naive -x \"0:2<umi1[6]>\" $test_dir/test.fq" e33034842bc6d691e4b84cb7943e76a2
checkcmdoutput "$splitcode --trim-only --pipe -q 10 --qtrim-3 --qtrim-naive --qtrim-pre -x \"0:2<umi1[6]>\" $test_dir/test.fq" cdf3d28488e3b6feb23efbe104281684
checkcmdoutput "$splitcode --trim-only --pipe -q 10 --qtrim-3 --qtrim-pre -x \"0:2<umi1[6]>\" $test_dir/test.fq" c883d391b1d4c20deac071f7f7af1905
checkcmdoutput "$splitcode --trim-only --pipe -q 10 --qtrim-3 --qtrim-pre -x \"0:2<umi1[3]>\" $test_dir/test.fq" 6b769e9893bd3e7f5a7cde91cf6420ed
checkcmdoutput "$splitcode --trim-only --pipe -q 10 --qtrim-3 --qtrim-pre -x \"0:2<umi1>0:-1\" -5 1 $test_dir/test.fq" 8626fc01d5b7db07ae63e3183b6d21e4
checkcmdoutput "$splitcode --trim-only --pipe -q 10 --qtrim-3 --qtrim-pre -5 4 -E ATCG $test_dir/test.fq" 4cdf8aa6f1767a22ad8c09e38fe2118b
checkcmdoutput "$splitcode --trim-only --pipe -q 10 --qtrim-3 --qtrim-pre -5 5 -E ATCG $test_dir/test.fq" 1b5a09bd343382ee78c9aa51245557c2
checkcmdoutput "$splitcode --trim-only --pipe -q 10 --qtrim-3 -w 5 $test_dir/test.fq" b0e38a846bb762767fcfa6ca002a7cab
checkcmdoutput "$splitcode --trim-only --pipe -q 10 --qtrim-3 -w 0:4 $test_dir/test.fq" ba907ffe8c0e14a57043a926f3ebce6b
checkcmdoutput "$splitcode --trim-only --pipe -q 10 --qtrim-3 --qtrim-pre -x \"0:2<umi1>0:-1\" -5 1 -3 3 $test_dir/test.fq" d90efb54e48c2c1c261cd4eb5e86cfdc
checkcmdoutput "$splitcode --trim-only --pipe -q 10 --qtrim-3 --qtrim-pre -x \"0:2<umi1>0:-1\" -5 1 -3 2 $test_dir/test.fq" 8626fc01d5b7db07ae63e3183b6d21e4
checkcmdoutput "$splitcode --trim-only --pipe -q 10 --qtrim-3 -x \"0:2<umi1>0:-1\" -5 1 -3 2 $test_dir/test.fq" 0d99695915b46e66e0dace6ce33fed44
checkcmdoutput "$splitcode --trim-only --pipe -q 10 --qtrim-3 -s $test_dir/test_summary.txt -x \"0:2<umi1>0:-1\" -5 1 -3 2 $test_dir/test.fq" 0d99695915b46e66e0dace6ce33fed44
cmdexec "wc $test_dir/test_summary.txt"
checkcmdoutput "$splitcode --trim-only --pipe -q 10 --qtrim-3 --qtrim-pre -5 5 -E ATCG $test_dir/test.fq" 1b5a09bd343382ee78c9aa51245557c2

# Adapter trimming tests

checkcmdoutput "$splitcode --trim-only -b CCAAA --partial5=3:0.35 --left=1 --pipe $test_dir/test.fq" b637fbabe71eb90bb9b3399a17eabef7
checkcmdoutput "$splitcode --trim-only -b CCAAA --partial5=3:0.34 --left=1 --pipe $test_dir/test.fq" b637fbabe71eb90bb9b3399a17eabef7
checkcmdoutput "$splitcode --trim-only -b CCAAA --partial5=3:0.33 --left=1 --pipe $test_dir/test.fq" cb52b79ed7469ca2ffe5739ec544b157
checkcmdoutput "$splitcode --trim-only -b CCAAA --partial5=4:0.34 --left=1 --pipe $test_dir/test.fq" cb52b79ed7469ca2ffe5739ec544b157
checkcmdoutput "$splitcode --trim-only -b CCAAA --partial5=2:0.34 --left=1 --pipe $test_dir/test.fq" c6eba12c36e53301f23a9823c2901f24
checkcmdoutput "$splitcode --trim-only -b CCAAA,CCGGAA --partial5=2:0.34, --partial3=,4 --left=1,0 --right=0,1 --pipe $test_dir/test.fq" 5d4541fb96da328d07ab9189216cf4a5
checkcmdoutput "$splitcode --trim-only -b CCGC -l 0:-4:0 --partial3=4:0.25 --right=1 --pipe $test_dir/test.fq" 11b55a195b5976331305569416db5bd4
checkcmdoutput "$splitcode --trim-only -b CCGG,CCGG -i a,b -l 0:-4:9,1:-4:10 --partial3=2,2 --right=1,1 -N 2 --pipe $test_dir/test.fq $test_dir/test.fq" 6807e3ba911fde8fb437f693d055c11f
checkcmdoutput "$splitcode --trim-only -b CCGC,GAAG -a ,{CCGC} -v {GAAG}, -l 0:-4:0, --partial3=4:0.25,3 --partial5=4,3 --right=1,0 --left=0,1 --pipe $test_dir/test.fq" 93f1726415edb410d5e733603bc4be11
checkcmdoutput "$splitcode --trim-only -b CCGC,GAAG -a ,{CCGC} -v {GAAG}, -l 0:-4:0, --partial3=4:0.25,3 --partial5=4,4 --right=1,0 --left=0,1 --pipe $test_dir/test.fq" cb52b79ed7469ca2ffe5739ec544b157
checkcmdoutput "$splitcode --trim-only -b CCGC,GAAG -a ,{CCGC} -l 0:-4:0, --partial3=4:0.25,3 --partial5=4,4 --right=1,0 --left=0,1 --pipe $test_dir/test.fq" 11b55a195b5976331305569416db5bd4

# More tests

echo "@read1
GATACTCCGGTTTTTTTTTTAAAAACCCCCGATACTCCGG
+
;I;<)(,%#\$KI;<)(,%#\$KI;<)(,%#\$\$\$\$\$\$\$\$\$\$\$
@read2
ATGCTTCCGGTTTTTTTTTTAAAAACCCCC
+
KI;<)(,%#\$KI;<)(,%#\$KI;<)(,%#\$
@read3
TAGAGAACACTTTTTTTTTTAAAAACCCCC
+
KI;<)(,%#\$KI;<)(,%#\$KI;<)(,%#\$
@read4
TTAATTCCGGTTTTTTTTTTAAAAACCCCC
+
KI;<)(,%#\$KI;<)(,%#\$KI;<)(,%#\$
@read5
TTGGTACTCCTTTTTTTTTTAAAAACCCCC
+
KI;<)(,%#\$KI;<)(,%#\$KI;<)(,%#\$
@read6
CCTAGAGAATTTTTTTTTTTAAAAACCCCC
+
KI;<)(,%#\$KI;<)(,%#\$KI;<)(,%#\$
@read7
GGTACTCCAATTTTTTTTTTAAAAACCCCC
+
KI;<)(,%#\$KI;<)(,%#\$KI;<)(,%#\$
@read8
TAGAGAACCTTTTTTTTTTTAAAAACCCCC
+
KI;<)(,%#\$KI;<)(,%#\$KI;<)(,%#\$" > $test_dir/test_1.fq

echo "@read1
GATACTCCGGTTTTTTTTTTAAAAACCCCC
+
KI;<)(,%#\$KI;<)(,%#\$KI;J;;;;;;
@read2
ATGCTTCCGGTTTTTTTTTTAAAAACCCCC
+
KI;<)(,%#\$KI;<)(,%#\$KI;<)(,%#\$
@read3
TAGAGAACACTTTTTTTTTTAAAAACCCCC
+
KI;<)(,%#\$KI;<)(,%#\$KI;<)(,%#\$
@read4
TTAATTCCGGTTTTTTTTTTAAAAACCCCC
+
KI;<)(,%#\$KI;<)(,%#\$KI;<)(,%#\$
@read5
TTGGTACTCCTTTTTTTTTTAAAAACCCCC
+
KI;<)(,%#\$KI;<)(,%#\$KI;<)(,%#\$
@read6
CCTAGAGAATTTTTTTTTTTAAAAACCCCC
+
KI;<)(,%#\$KI;<)(,%#\$KI;<)(,%#\$
@read7
GGTACTCCAATTTTTTTTTTAAAAACCCCC
+
KI;<)(,%#\$KI;<)(,%#\$KI;<)(,%#\$
@read8
TAGAGAACCTTTTTTTTTTTAAAAACCCCC
+
KI;<)(,%#\$KI;<)(,%#\$KI;<)(,%#\$" > $test_dir/test_2.fq

# Testing sam tags (and some other checks)

checkcmdoutput "$splitcode --pipe --no-outb  --x-names --seq-names --com-names --nFastqs=2 --keep-grp=<(echo "t,b") -m /dev/null -g t,t,b,b,b -b GATACT,TAGAGAA,CCGG,AACC,AAAAACC -d 1,1,1,1,1 -l 0:0:10,0:0:10,0,0,0 -v ,,{{t}}14-14,{{t}}14-14,{{t}}14-14 $test_dir/test_1.fq $test_dir/test_2.fq" 53a53635d0ad34937feaca1f2aa8ae5d
checkcmdoutput "$splitcode --pipe --no-outb  --x-names --seq-names --com-names --nFastqs=2 --keep-grp=<(echo "t,b") -m /dev/null -g t,t,b,b,b -b GATACT,TAGAGAA,CCGG,AACC,AAAAACC -d 1,1,1,1,1 -l 0:0:10,0:0:10,0,0,0 -v ,,{{t}}14-14,{{t}}14-14,{{t}}14-14 -M AB:G:,XX:X:,FS:D:,JKL $test_dir/test_1.fq $test_dir/test_2.fq" fd979c64c21eb530e3fb8eb964efac25
checkcmdoutput "$splitcode --pipe --no-outb  --x-names --seq-names --com-names --nFastqs=2 -m /dev/null -g t,t,b,b,b -b GATACT,TAGAGAA,CCGG,AACC,AAAAACC -d 0,0,1,1,1 -l 0:0:10,0:0:7,0,0,0 -v ,,{{t}}14-14,{{t}}15-16,{{t}}15-16 -M AB:G:,XX:X:,FS:D:,JKL $test_dir/test_1.fq $test_dir/test_2.fq" 6620fed5ace9d1a1a85a51065d5b0c11

# Testing subs

checkcmdoutput "$splitcode --pipe --no-outb  --x-names --seq-names --com-names --nFastqs=2 --keep-grp=<(echo "t,b") -m /dev/null -g t,t,b,b,b -b GATACT,TAGAGAA,CCGG,AACC,AAAAACC -d 1,1,1,1,1 -l 0:0:10,0:0:10,0,0,0 -v ,,{{t}}14-14,{{t}}14-14,{{t}}14-14 -M AB:G:,XX:X:,FS:D:,JKL --subs=.,.,.,.,. $test_dir/test_1.fq $test_dir/test_2.fq" 8e9a69d297d492cbfe054004666efa44
checkcmdoutput "$splitcode --pipe --no-outb  --x-names --seq-names --com-names --nFastqs=2 --keep-grp=<(echo "t,b") -m /dev/null -g t,t,b,b,b -b GATACT,TAGAGAA,CCGG,AACC,AAAAACC -d 1,1,1,1,1 -l 0:0:10,0:0:10,0,0,0 -M AB:G:,XX:X:,FS:D:,JKL --subs=.,.,.,.,. $test_dir/test_1.fq $test_dir/test_2.fq" 078e57e27856ccb38ca751a9a832bfc3
checkcmdoutput "$splitcode --pipe --no-outb  --x-names --seq-names --com-names --nFastqs=2 -m /dev/null -g t,t,b,b,b -b GATACT,TAGAGAA,CCGG,AACC,AAAAACC -d 1,1,1,1,1 -l 0:0:10,0:0:10,0,0,0 -M AB:G:,XX:X:,FS:D:,JKL --subs=.,.,.,.,. $test_dir/test_1.fq $test_dir/test_2.fq" 05940a828e3dbfbd8a0c44f9b375abb6

checkcmdoutput "$splitcode --trim-only --pipe --x-names --seq-names --nFastqs=2 --keep-grp=<(echo "t,b") -g t,t,b,b,b -b GATACT,TAGAGAA,CCGG,AACC,AAAAACC -d 1,1,1,1,1 -l 0:0:10,0:0:10,0,0,0 -v ,,{{t}}14-14,{{t}}14-14,{{t}}14-14 --subs=TT,GG,GG,GG,GG $test_dir/test_1.fq $test_dir/test_2.fq" 47191e11ace788948cea359400a66b55
checkcmdoutput "$splitcode --trim-only --pipe --x-names --seq-names --nFastqs=2 -g t,t,b,b,b -b GATACT,TAGAGAA,CCGG,AACC,AAAAACC -d 1,1,1,1,1 -l 0:0:10,0:0:10,0,0,0 --subs=TT,GG,GG,GG,GG $test_dir/test_1.fq $test_dir/test_2.fq" 1ef6f5045ebd99f29afae768ecbf79f3
checkcmdoutput "$splitcode --trim-only --pipe --x-names --seq-names --nFastqs=2 -g t,t,b,b,b -b GATACT,TAGAGAA,CCGG,AACC,AAAAACC -d 1,1,1,1,1 -l 0:0:10,0:0:10,0,0,0 --subs=TT,-,GG,GG,GG $test_dir/test_1.fq $test_dir/test_2.fq" c2abb7ce085324c0bbb89ae9f01dd85d
checkcmdoutput "$splitcode --trim-only --pipe --x-names --seq-names --nFastqs=2 -b GATACT,TAGAGAA,CCGG,AACC,AAAAACC,T:5-11 -d 1,1,1,1,1,0 -l 0:0:10,0:0:10,0,0,0,1 --subs=TT,GG,GG,GG,GG,- $test_dir/test_1.fq $test_dir/test_2.fq" 6b3e2fb7655a14b08cddf7effd5f27e8
checkcmdoutput "$splitcode --trim-only --pipe --x-names --seq-names --nFastqs=1 -b T:5-6,C:3-100,A:1-5 --subs=.,-,- $test_dir/test_1.fq" 8cc70e8d1a5bf62d602c38558ecbd725
checkcmdoutput "$splitcode --trim-only --pipe --x-names --seq-names --nFastqs=1 -b ACAC,C:5-100,A:1-5 --subs=-,-,- $test_dir/test_1.fq" 3bb7150085b56c492f58d6112e234bd5
checkcmdoutput "$splitcode --trim-only --pipe --x-names --seq-names --nFastqs=2 -b GATACT,TAGAGAA,CCGG,AACC,AAAAACC,AC:1-2 -d 1,1,1,1,1,0 -l 0:0:10,0:0:10,0,0,0,1 --subs=TT,GG,GG,GG,GG,- $test_dir/test_1.fq $test_dir/test_2.fq" b6d90e8d48b3bd9c6939c8b095b18804
checkcmdoutput "$splitcode --trim-only --pipe --x-names --seq-names --nFastqs=2 -b GATACT,TAGAGAA,CCGG,AACC,AAAAACC,AC:2 -d 1,1,1,1,1,0 -l 0:0:10,0:0:10,0,0,0,1 --subs=TT,GG,GG,GG,GG,- $test_dir/test_1.fq $test_dir/test_2.fq" 83d5094e21e3044c9e8ea03e13aac2cd
checkcmdoutput "$splitcode --trim-only --pipe --x-names --seq-names --nFastqs=2 -b GATACT,TAGAGAA,CCGG,AACC,AAAAACC,AC:2 -d 1,1,1,1,1,0 -l 0:0:10,0:0:10,0,0,0,1 --subs=TT,GG,GG,GG,GG,NNNNNNNN $test_dir/test_1.fq $test_dir/test_2.fq" d5f194caf779d8a1fd0aa59bb9fe8eae

checkcmdoutput "$splitcode --trim-only --pipe --x-names --seq-names --nFastqs=2 -g t,t,b,b,b -b GATACT,TAGAGAA,CCGG,AACC,AAAAACC -d 1,1,1,1,1 -v ,,{{t}}14-14,{{t}}14-14,{{t}}14-14 --subs=TT,GG,GG,GG,GG $test_dir/test_1.fq $test_dir/test_2.fq" 4d8a551a84d4d10b7deac330227616cc
checkcmdoutput "$splitcode --trim-only --pipe --x-names --seq-names --nFastqs=2 -g t,t,b,b,b -b GATACT,TAGAGAA,CCGG,AACC,AAAAACC -d 1,1,1,1,1 -v ,,{{t}}14-14,{{t}}14-14,{{t}}14-14 --subs=TTTTTTTTTTTTTTT,GGGGGGGGGGGGGGG,GG,GG,GG -q 30 --qtrim-5 -M AB:G:,XX:X:,FS:D:,JKL -n 2 $test_dir/test_1.fq $test_dir/test_2.fq" 622e99643d60536a6ee99f9d968ac20f
checkcmdoutput "$splitcode --trim-only --pipe --x-names --seq-names --nFastqs=2 -g t,t,b,b,b -b GATACT,TAGAGAA,CCGG,AACC,AAAAACC -d 1,1,1,1,1 -v ,,{{t}}14-14,{{t}}14-14,{{t}}14-14 --subs=TTTTTTTTTTTTTTT,GGGGGGGGGGGGGGG,GG,GG,GG -q 30 --qtrim-3 --qtrim-5 -n 1 $test_dir/test_1.fq $test_dir/test_2.fq" 91071e8bbcecbad10599a7d8ebd3815c
checkcmdoutput "$splitcode --trim-only --pipe --x-names --seq-names --nFastqs=2 -g t,t,b,b,b -b GATACT,TAGAGAA,CCGG,AACC,AAAAACC -d 1,1,1,1,1 -v ,,{{t}}14-14,{{t}}14-14,{{t}}14-14 --subs=TTTTTTTTTTTTTTT,GGGGGGGGGGGGGGG,GG,GG,GG -q 30 --qtrim-5 --qtrim-pre -M AB:G:,XX:X:,FS:D:,JKL -n 2 $test_dir/test_1.fq $test_dir/test_2.fq" 901358fbfe41d20a9b4dd99fbac22866
checkcmdoutput "$splitcode --trim-only --pipe --x-names --seq-names --nFastqs=2 -g t,t,b,b,b -b GATACT,TAGAGAA,CCGG,AACC,AAAAACC -d 1,1,1,1,1 -v ,,{{t}}14-14,{{t}}14-14,{{t}}14-14 --subs=TTTTTTTTTTTTTTT,GGGGGGGGGGGGGGG,GG,GG,GG -q 30 --qtrim-5 --qtrim-pre -3 2,1 -M AB:G:,XX:X:,FS:D:,JKL -n 2 $test_dir/test_1.fq $test_dir/test_2.fq" e1dbf41e8fbc3b9919834a241feb9bb8
checkcmdoutput "$splitcode --trim-only --pipe --x-names --seq-names --nFastqs=2 -g t,t,b,b,b -b GATACT,TAGAGAA,CCGG,AACC,AAAAACC -d 1,1,1,1,1 -v ,,{{t}}14-14,{{t}}14-14,{{t}}14-14 --subs=TTTTTTTTTTTTTTT,GGGGGGGGGGGGGGG,GG,GG,GG -5 3,2 -3 2,1 -M AB:G:,XX:X:,FS:D:,JKL -n 2 $test_dir/test_1.fq $test_dir/test_2.fq" 26c827664525506605c10cbfe6ce8566
checkcmdoutput "$splitcode --trim-only --pipe --x-names --seq-names --nFastqs=2 -g t,t,b,b,b -b GATACT,TAGAGAA,CCGG,AACC,AAAAACC -d 1,1,1,1,1 -v ,,{{t}}14-14,{{t}}14-14,{{t}}14-14 --subs=TTTTTTTTTTTTTTT,GGGGGGGGGGGGGGG,GG,GG,GG -n 3 --left=1:1,1,0,0,0 $test_dir/test_1.fq $test_dir/test_2.fq" 867987cdb5a7fec25ae37f6d2742aafd
checkcmdoutput "$splitcode --trim-only --pipe --x-names --seq-names --nFastqs=2 -g t,t,b,b,b -b GATACT,TAGAGAA,CCGG,AACC,AAAAACC -d 1,1,1,1,1 -v ,,{{t}}14-14,{{t}}14-14,{{t}}14-14 --subs=TTTTTTTTTTTTTTT,GGGGGGGGGGGGGGG,GG,GG,GG -n 3 --left=1:1,1,1,1,1 $test_dir/test_1.fq $test_dir/test_2.fq" dd3e052e2383110ae6fbfb256496efe2
checkcmdoutput "$splitcode --trim-only --pipe --x-names --seq-names --nFastqs=2 -g t,t,b,b,b -b GATACT,TAGAGAA,CCGG,AACC,AAAAACC -d 1,1,1,1,1 -v ,,{{t}}14-14,{{t}}14-14,{{t}}14-14 --subs=TTTTTTTTTTTTTTT,GGGGGGGGGGGGGGG,GG,GG,GG -n 3 --left=1:1,1,1,1,1 -q 30 --qtrim-5 -w 1,1 $test_dir/test_1.fq $test_dir/test_2.fq" 7e6de423bcff3285977a599530e0e1fd


checkcmdoutput "$splitcode -m /dev/null --pipe --x-names --seq-names --nFastqs=2 -g t,t,b,b,b -b GATACT,TAGAGAA,CCGG,AACC,AAAAACC -d 1,1,1,1,1 -v ,,{{t}}14-14,{{t}}14-14,{{t}}14-14 --subs=TTTTTTTTTTTTTTT,GGGGGGGGGGGGGGG,GG,GG,GG -q 30 --qtrim-5 --qtrim-pre -3 2,1 -M AB:G:,XX:X:,FS:D:,JKL -n 2 --extract=\"{GATACT}<umi[2]>\" --no-outb --no-chain --no-x-out -l 0,0,0,0,0 $test_dir/test_1.fq $test_dir/test_2.fq" 0f55956a9552e59fa5dfde6cf9d758ad
checkcmdoutput "$splitcode -m /dev/null --pipe --x-names --seq-names --nFastqs=2 -g t,t,b,b,b -b GATACT,TAGAGAA,CCGG,AACC,AAAAACC -d 1,1,1,1,1 -v ,,{{t}}14-14,{{t}}14-14,{{t}}14-14 --subs=TTTTTTTTTTTTTTT,GGGGGGGGGGGGGGG,GG,GG,GG -q 30 --qtrim-5 --qtrim-pre -3 2,1 -M AB:G:,XX:X:,FS:D:,JKL -n 2 --extract=\"{GATACT}<umi[3]>\" --no-outb --no-chain --no-x-out -l 0,0,0,0,0 $test_dir/test_1.fq $test_dir/test_2.fq" 022da965571c430ad0bd98455cb21820

checkcmdoutput "$splitcode -m $test_dir/testmapping.txt --pipe --x-names --seq-names --nFastqs=2 -g t,t,b,b,b -b GATACT,TAGAGAA,CCGG,AACC,AAAAACC -d 1,1,1,1,1 -v ,,{{t}}14-14,{{t}}14-14,{{t}}14-14 --subs=TTTTTTTTTTTTTTT,GGGGGGGGGGGGGGG,GG,GG,GG -n 3 --left=1:1,1,1,1,1 -q 30 --qtrim-5 --mod-names -w 1,1 $test_dir/test_1.fq $test_dir/test_2.fq" d41d8cd98f00b204e9800998ecf8427e
checkcmdoutput "cat $test_dir/testmapping.txt" adee4c7edfa0b9bb26d6ebb4d351582d

# Additional advanced UMI extraction tests (special cases)

checkcmdoutput "$splitcode -t 1 -N 2 -c $test_dir/splitcode_example_config.txt -y <(echo \"DPM,Y,ODD,EVEN,ODD\") --x-names -m /dev/null --pipe -x \"<umi1{{EVEN}}>,<umi2{{ODD}}>\" $test_dir/A_1.fastq.gz $test_dir/A_2.fastq.gz|head -20" 4dc03f395f130585d8fb9b6d97ec8efa
checkcmdoutput "$splitcode -t 1 -N 2 -c $test_dir/splitcode_example_config.txt -y <(echo \"DPM,Y,ODD,EVEN,ODD\") --x-names -m /dev/null --pipe -x \"<umi1{{EVEN}}>,<umi2{{ODD}}>\" $test_dir/A_1.fastq.gz $test_dir/A_2.fastq.gz|head -560|tail -20" 42243836798d98d3dd1cb9b054029929
checkcmdoutput "$splitcode -t 1 --no-chain -N 2 -c $test_dir/splitcode_example_config.txt -y <(echo \"DPM,Y,ODD,EVEN,ODD\") --x-names -m /dev/null --pipe -x \"<umi1{{EVEN}}>,<umi2{{ODD}}>\" $test_dir/A_1.fastq.gz $test_dir/A_2.fastq.gz|head -560|tail -20" eaba520366b8572afac57539d69e716e
checkcmdoutput "$splitcode -t 1 -N 2 -c $test_dir/splitcode_example_config.txt -y <(echo \"DPM,Y,ODD,EVEN,ODD\") --x-names -m /dev/null --pipe -x \"1:10<umi1>1:20,1:10<umi2[8-12]>1:20,<umi1{{EVEN}}>,<umi2{{ODD}}>,<umi3{{ODD}}>,<umi4{Odd2Bo10}>\" $test_dir/A_1.fastq.gz $test_dir/A_2.fastq.gz|head -20" d6d98eb83e4e60e3a9f1ace688db065a
checkcmdoutput "$splitcode --trim-only -t 1 -i bc1 -b AAAAAGCCCC -d 1 --x-names --pipe -x \"<~umi1{bc1}>\" $test_dir/test_1.fq" 51e0dc9135b51a3ec7399882913ccb99
cmdexec "$splitcode --trim-only -t 1 -i bc1 -b AAAAAGCCCC -d 1 -o /dev/null -x \"<~$test_dir/umi1{bc1}>\" $test_dir/test_1.fq"
checkcmdoutput "cat $test_dir/umi1.fastq" 18b0a455896b74aab83a6b0ff721de05
checkcmdoutput "$splitcode --trim-only -t 1 -i bc1 -b AAAAAGCCCC -d 1 --x-names --pipe --subs=ATC -x \"<~umi1{#bc1}>\" $test_dir/test_1.fq" 56a0fb91eb32570b511e0c9525263a3a
checkcmdoutput "$splitcode --trim-only -t 1 -i bc1 -b AAAAAGCCCC -d 1 --x-names --pipe -x \"<~umi1{@bc1}>\" $test_dir/test_1.fq" 8b434ce88a9d5bb6af4743b0cccbf209
checkcmdoutput "$splitcode --trim-only -t 1 -i bc1 -b AAAAAGCCCC -d 1 --x-names --pipe -x \"<~umi1{@bc1}>,<umi2{bc1}>,<umi1[3]>{bc1}\" $test_dir/test_1.fq" 813e94a0fa68aa1a9ea2c3a74b7daed6
checkcmdoutput "$splitcode --trim-only -t 1 -i bc1 -b AAAAAGCCCC -d 1 --x-names --pipe -x \"<~umi1{@bc1}>,<umi2{bc1}>,<umi1[3]>{bc1},<umi3[3]>{bc1}\" --sam-tags=\"CB:Z:,RX:Z:/RY:Z:\" $test_dir/test_1.fq" e675ac734c9859c0a03a5bbc967a0c6f
checkcmdoutput "$splitcode -t 1 -N 2 -c $test_dir/splitcode_example_config.txt -y <(echo \"DPM,Y,ODD,EVEN,ODD\") --x-names -m /dev/null --pipe -x \"<umi1{{EVEN}}>,<~umi3{*}>,<umi2{{ODD}}>\" $test_dir/A_1.fastq.gz $test_dir/A_2.fastq.gz|head -48" b340f4e12eeaddffa88ffc83c6964afa
checkcmdoutput "$splitcode -t 1 -N 2 -c $test_dir/splitcode_example_config.txt -y <(echo \"DPM,Y,ODD,EVEN,ODD\") --x-names -m /dev/null --pipe -x \"<umi1{{EVEN}}>,<umi2{*}>,<umi2{{ODD}}>\" $test_dir/A_1.fastq.gz $test_dir/A_2.fastq.gz|head -40" 130d75560550456801354c9b31ba4a27
checkcmdoutput "$splitcode --trim-only -t 1 -i bc1 -b AAAAACCCCCG --x-names --pipe -x \"0:0<umi1[5]>,<~umi2{@bc1}>,<~umi3{*}>\" --no-x-out --sam-tags=\"CB:Z:,RX:Z:/RY:Z:/RZ:Z:,BI:I:\" $test_dir/test_1.fq" b7cb8dd7e44b34f07f66ea3b3a9e7504
checkcmdoutput "$splitcode --trim-only -t 1 -i bc0,bc1 -b TACT,AAAAACCCCCG -g a,a --x-names --pipe -x \"0:0<umi1[5]>,<~umi2[4]{{a}}>,<~umi3[15-25]{*}>\" --no-x-out --sam-tags=\"CB:Z:,RX:Z:/RY:Z:/RZ:Z:,BI:I:\" $test_dir/test_1.fq" f06262837c002c9e37e28535d2a617ae
checkcmdoutput "$splitcode --trim-only -t 1 -i bc0,bc1 -b TACT,AAAAACCCCCG -g a,a --subs=AA,TT --x-names --pipe -x \"0:0<umi1[5]>,<~umi2{{a}}>,<~umi3{#*}>\" --no-x-out --sam-tags=\"CB:Z:,RX:Z:/RY:Z:/RZ:Z:,BI:I:\" $test_dir/test_1.fq" 4cb0d9c554a74b9075060ce6c2e705e9
checkcmdoutput "$splitcode --trim-only -t 1 -i bc0,bc1 -b TACT,AAAAACCCCCGC -d 0,1 -g a,a --subs=AA,TT --x-names --pipe -x \"0:0<umi1[5]>,<~umi2{{a}}>,<~umi3{@*}>\" --no-x-out --sam-tags=\"CB:Z:,RX:Z:/RY:Z:/RZ:Z:,BI:I:\" -N 2 $test_dir/test_1.fq $test_dir/test_1.fq" 3cf1eec93700913db553718a7440df03

# Testing --sub-assign

checkcmdoutput "$splitcode -t 1 --pipe -N 2 -c $test_dir/splitcode_example_config.txt -y <(echo "DPM,Y,ODD,EVEN,ODD") --com-names -m /dev/null --sub-assign=5 $test_dir/A_1.fastq.gz $test_dir/A_2.fastq.gz" d9bd0b3183cfb50256d5473af2f6476e
checkcmdoutput "$splitcode -t 1 --pipe -N 2 -c $test_dir/splitcode_example_config.txt -y <(echo "DPM,Y,ODD,EVEN,ODD") --com-names -m /dev/null --sub-assign=4,5,6 $test_dir/A_1.fastq.gz $test_dir/A_2.fastq.gz" c2b234369f465fc982687bb71c3750c4
checkcmdoutput "$splitcode -t 1 --pipe -N 2 -c $test_dir/splitcode_example_config.txt -y <(echo "DPM,Y,ODD,EVEN,ODD") --com-names -m /dev/null --sub-assign=0 $test_dir/A_1.fastq.gz $test_dir/A_2.fastq.gz" b8a61b0af84661c4da6fb4c7151a2210
checkcmdoutput "$splitcode -t 1 --pipe -N 2 -c $test_dir/splitcode_example_config.txt -y <(echo "DPM,Y,ODD,EVEN,ODD") --com-names -m /dev/null --sub-assign=1,2 $test_dir/A_1.fastq.gz $test_dir/A_2.fastq.gz" c4f74a40c002f329d1871f4611aadd92
checkcmdoutput "$splitcode -t 1 --pipe -N 2 -c $test_dir/splitcode_example_config.txt -y <(echo "DPM,Y,ODD,EVEN,ODD") --com-names -m /dev/null --sub-assign=0,1,2 $test_dir/A_1.fastq.gz $test_dir/A_2.fastq.gz" c90450075bda6beab6656095fb344991
checkcmdoutput "$splitcode -t 1 --pipe -N 2 -c $test_dir/splitcode_example_config.txt -y <(echo \"DPM,Y,ODD,EVEN,ODD\") --com-names -m /dev/null --sub-assign=0,1,2,3,4 --sam-tags=\"CB:Z:,RX:Z:,BI:i:,AA:i:\" $test_dir/A_1.fastq.gz $test_dir/A_2.fastq.gz" 9af079a5eef6bc1f35b438a39b48b618

# Testing piping and --keep-com

checkcmdoutput "$splitcode --assign -t 1 --pipe -N 2 -c $test_dir/splitcode_example_config.txt -y <(echo \"DPM,Y,ODD,EVEN,ODD\") --mod-names --com-names --no-outb -m /dev/null $test_dir/A_1.fastq.gz $test_dir/A_2.fastq.gz|head -96|$splitcode --assign -t 1 --pipe -N 2 -c $test_dir/splitcode_example_config.txt -y <(echo \"DPM,Y,ODD,EVEN,ODD\") --com-names --keep-com --sub-assign=0 -m /dev/null --inleaved -" 0e7e7a8f68481e9262f680068c14f971

# Testing remultiplexing

echo "AAA $test_dir/A_1.fastq.gz $test_dir/A_2.fastq.gz
BB $test_dir/B_1.fastq.gz $test_dir/B_2.fastq.gz
CC $test_dir/test.fq $test_dir/test.fq
DD $test_dir/test_1.fq $test_dir/test_2.fq" > $test_dir/batch.txt

checkcmdoutput "$splitcode --trim-only --remultiplex -t 1 --pipe --com-names $test_dir/batch.txt" cf242391e572712067750d6b4169918d
cmdexec "$splitcode --trim-only --prefix=G -c $test_dir/splitcode_example_config.txt -m $test_dir/ttmap.txt.gz -O $test_dir/testO.fq.gz -u \"$test_dir/ttu1.fq.gz,$test_dir/ttu2.fq.gz\" -x \"0:0<$test_dir/umix[2]>\" -o \"$test_dir/tt1.fq.gz,$test_dir/tt2.fq.gz\" --remultiplex -t 1 -y <(echo \"DPM,Y,ODD,EVEN,ODD\") --com-names --bc-names $test_dir/batch.txt"
checkcmdoutput "zcat < $test_dir/ttmap.txt.gz" 46ece172aa60278e06f68e38c41a0497
checkcmdoutput "zcat < $test_dir/testO.fq.gz" 643001abd6b528f5c629b17977bd1f14
checkcmdoutput "zcat < $test_dir/ttu2.fq.gz" 51611c4170e95d7787182e54666af4a8
checkcmdoutput "zcat < $test_dir/umix.fastq.gz" 56654b92926b5edf2f08e9b715c41167
checkcmdoutput "zcat < $test_dir/tt1.fq.gz" 668abd8a81136052c5ea1be491e23b9a

# Testing boundaries

echo "@read1
ATGAGT
+
FFFFFF
@read2
ACACTT
+
FFFFFF" > $test_dir/test_bound.fq

checkcmdoutput "$splitcode --trim-only -b GAG -p --mod-names -l "0:-40:0" -R 1 -x \"0:1<umi[4]>\" $test_dir/test_bound.fq" 4e32e75c064c7e4c7b2edf277ee839e1

# Some additional tests

checkcmdoutput "$splitcode --trim-only -b AT,TT --loc-names -d 1 -p --mod-names $test_dir/test_bound.fq" 432d31a38782e90fc8129c0d72e92d16
checkcmdoutput "$splitcode --trim-only -b AT,TT,GGGGG -i a,a,b --loc-names -d 1 -p --mod-names $test_dir/test_bound.fq" ca97fa58494f66bd3395ac329e5a95cf
checkcmdoutput "$splitcode --trim-only -b AT,TT,GGGGG -i a,a,a --loc-names -d 1 -p --mod-names $test_dir/test_bound.fq" ca97fa58494f66bd3395ac329e5a95cf

# Testing --min-delta

checkcmdoutput "$splitcode --trim-only -b AGG,ACC -d 2 --min-delta=0 -p --mod-names $test_dir/test_bound.fq" 08a1669f62d83b7726d562be030f44b4

checkcmdoutput "$splitcode --trim-only -b AGG,ACC -d 2 --min-delta=1 -p --mod-names $test_dir/test_bound.fq" af972eebe3b6e9366dacda3c0902af82

checkcmdoutput "$splitcode --trim-only -b AGG,ACC -d 2 --min-delta=2 -p --mod-names $test_dir/test_bound.fq" af972eebe3b6e9366dacda3c0902af82

checkcmdoutput "$splitcode --trim-only -b ATGAGA,GGGGGG,ATGACC -d 2 --min-delta=0 -p --mod-names $test_dir/test_bound.fq" decb5f9ccbd1300e9b89a78b52fc9228

checkcmdoutput "$splitcode --trim-only -b ATGAGA,GGGGGG,ATGACC -d 2 --min-delta=1 -p --mod-names $test_dir/test_bound.fq" 40e7d95afa0ffa0fc483fc87f0c8a998

checkcmdoutput "$splitcode --trim-only -b ATGAGA,GGGGGG,ATGACC -d 2 --min-delta=2 -p --mod-names $test_dir/test_bound.fq" 40e7d95afa0ffa0fc483fc87f0c8a998

checkcmdoutput "$splitcode --trim-only -b TGAGC,TGAGA,TGAGG -d 1 --min-delta=0 --mod-names -p $test_dir/test_bound.fq" 40e7d95afa0ffa0fc483fc87f0c8a998

checkcmdoutput "$splitcode --trim-only -b TGAGC,TGAGA,TGAGG -d 1 --min-delta=1 --mod-names -p $test_dir/test_bound.fq" 40e7d95afa0ffa0fc483fc87f0c8a998

checkcmdoutput "$splitcode --trim-only -b TGAGC,TGAGA,TGAGT -d 1 --min-delta=0 --mod-names -p $test_dir/test_bound.fq" 588c7374ce94f3e07d94b5630a9c3f4e

checkcmdoutput "$splitcode --trim-only -b TGAGC,TGAGA,TGAGT -d 1 --min-delta=1 --mod-names -p $test_dir/test_bound.fq" 40e7d95afa0ffa0fc483fc87f0c8a998



# Shorten assigned final barcode 
                
checkcmdoutput "$splitcode -m /dev/null -b AT,TT -d 1 -p --mod-names --bclen=15 --bc-names $test_dir/test_bound.fq" 9e7b7980f24349c5c334e1fc02950af8

# Try them with FASTA file

echo ">read1
ATGAGT
>read2
ACACTT" > $test_dir/test_bound.fa

checkcmdoutput "$splitcode --trim-only -b AT,TT --loc-names -d 1 -p --mod-names $test_dir/test_bound.fa" ec7150e2b71fcf683ed80a3013e7c203
checkcmdoutput "$splitcode --trim-only -b AT,TT,GGGGG -i a,a,b --loc-names -d 1 -p --mod-names $test_dir/test_bound.fa" 8d701ebb09497ab0a3031e22b3ebf8a0
checkcmdoutput "$splitcode --trim-only --out-fasta -b AT,TT --loc-names -d 1 -p --mod-names $test_dir/test_bound.fa" 3458061c02f7a942a1b4eddf47d5bd7a
checkcmdoutput "$splitcode --trim-only --out-fasta -b AT,TT,GGGGG -i a,a,b --loc-names -d 1 -p --mod-names $test_dir/test_bound.fa" 5295712caa2eab75f1464c3a6da88be8

# Allow multiple extracts in config file

echo "@extract 0:0<x[2]>,0:2<x[2]>" > $test_dir/config_mult_extracts.txt
echo "@trim-5 1" >> $test_dir/config_mult_extracts.txt
echo "@extract 0:3<x[1]>" >> $test_dir/config_mult_extracts.txt

checkcmdoutput "$splitcode --trim-only -c $test_dir/config_mult_extracts.txt --x-only -p $test_dir/test_bound.fq" aec8f4f4691db99827b80a0f21140848

echo "@no-chain" >> $test_dir/config_mult_extracts.txt
checkcmdoutput "$splitcode --trim-only -c $test_dir/config_mult_extracts.txt --x-only -p $test_dir/test_bound.fq" fad7c4fa3cb392c57255d2a8139d4ba0


echo "@extract 0:0<y[1]>" >> $test_dir/config_mult_extracts.txt
echo "@extract 0:3<y[1]>" >> $test_dir/config_mult_extracts.txt
echo "@no-chain y" >> $test_dir/config_mult_extracts.txt

checkcmdoutput "$splitcode --trim-only -c $test_dir/config_mult_extracts.txt --x-only -p $test_dir/test_bound.fq" 3bbf54b385da42b516117a9b1d0bc25d

echo "@no-chain x,y" >> $test_dir/config_mult_extracts.txt

checkcmdoutput "$splitcode --trim-only -c $test_dir/config_mult_extracts.txt --x-only -p $test_dir/test_bound.fq" 676235664622f341a6fcb6b2035642f2

# Check variable length tags across different locations

checkcmdoutput "$splitcode --trim-only --mod-names -b TCCT,TCCC,TCCA,TCC,TCCGG,AAA --previous=,,,{AAA},, -l 0:0:5,0:0:5,0:0:5,0,0:5:10,0:5:10 -p $test_dir/test.fq" a2df9f98d58b47aa1beeebabc835b3b3

# Check --unmask

echo ">read1
ATGAGT
>read2
ACACTT" > $test_dir/test_unmasked.1.fa

echo ">read1
ANNNNT
>read2
NCACTN" > $test_dir/test_unmasked.2.fa

checkcmdoutput "$splitcode --unmask $test_dir/test_unmasked.1.fa $test_dir/test_unmasked.2.fa \"\"" b98fba82fc4a92c49d916ab97a8463bc

# Check @remove

echo "@remove:" > $test_dir/config_with_remove.txt
echo "tag1" >> $test_dir/config_with_remove.txt
echo "" >> $test_dir/config_with_remove.txt
echo "id tag" >> $test_dir/config_with_remove.txt
echo "tag1 AAA" >> $test_dir/config_with_remove.txt
echo "@read1" > $test_dir/remove_input.fastq
echo "AAA" >> $test_dir/remove_input.fastq
echo "+" >> $test_dir/remove_input.fastq
echo ";;;" >> $test_dir/remove_input.fastq

checkcmdoutput "$splitcode --trim-only -c $test_dir/config_with_remove.txt --mod-names --out-fasta --nFastqs=1 --pipe $test_dir/remove_input.fastq" d41d8cd98f00b204e9800998ecf8427e

# Check barcode-encode in config file

echo "@assign" > $test_dir/bcencode.txt
echo "@barcode-encode X,X,X,X,X,X,X,X,X,X,X,X,X,X,X,X" >> $test_dir/bcencode.txt
echo "" >> $test_dir/bcencode.txt
echo "tags groups" >> $test_dir/bcencode.txt
echo "A X" >> $test_dir/bcencode.txt
echo "T X" >> $test_dir/bcencode.txt
echo "C X" >> $test_dir/bcencode.txt
echo "G X" >> $test_dir/bcencode.txt

echo "@read1" > $test_dir/bcencode.fastq
echo "AAAAAAAAAAAAAAAA" >> $test_dir/bcencode.fastq
echo "+" >> $test_dir/bcencode.fastq
echo "AAAAAAAAAAAAAAAA" >> $test_dir/bcencode.fastq

checkcmdoutput "$splitcode -m /dev/null --bclen=30 -c $test_dir/bcencode.txt --mod-names --out-fasta --nFastqs=1 --pipe $test_dir/bcencode.fastq" 6c300074e39d272b4fee6a75fae3e8b4

# Check extract ^...^ and ^^...^^

checkcmdoutput "$splitcode --trim-only -x \"0:1<^AAA^umi[3]>\" -p $test_dir/test_bound.fa" 5827f14f608ed69bd38b12602c8fdeea
checkcmdoutput "$splitcode --trim-only -x \"0:1<~^^CGC^^umi[3]>\" -p $test_dir/test_bound.fa" 20b95054c5a91cfeb8585c550570bbe1

# Check early termination (fixes a bug when maxFindsG is not set for certain groups)

echo "@rNYBot1_StgOdd2Bo1Part
TATTATGGTCCCCCCCCTTCGTGGAATCTAGCAAAAAATGACTTGAGCTGGGGGGGGGG
+
TATTATGGTCCCCCCCCTTCGTGGAATCTAGCTGACTTGTGACTAGCTDGGGGGGGGGG" > $test_dir/test_term.fq

checkcmdoutput "$splitcode --assign  -m /dev/null --mod-names -g PART,RPM,RPM,ODD,Y -i Part,RTBC2,RTBC4,Odd2Bo1,NYBot1_Stg -b TGACTTG,TTTTTTT,GGGGGGG,TTCGTGGAATCTAGC,TATTATGGT --maxFindsG=Y:1 -p $test_dir/test_term.fq" 79ad1ca9184b274e1396b5f2220cbade


