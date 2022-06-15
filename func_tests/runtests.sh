#!/bin/bash

### Setup ###

splitcode="./src/splitcode"
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

cmdexec "$splitcode --version"

cmdexec "$splitcode -t 1 -o "$test_dir/A_out_1.fastq.gz,$test_dir/A_out_2.fastq.gz" -O $test_dir/A_out_barcodes.fastq.gz -N 2 -c $test_dir/splitcode_example_config.txt -y <(echo "DPM,Y,ODD,EVEN,ODD") --mod-names --gzip -m $test_dir/A_out_mapping.txt.gz $test_dir/A_1.fastq.gz $test_dir/A_2.fastq.gz"
checkcmdoutput "zcat < $test_dir/A_out_1.fastq.gz" 05eb0a03078777fb36ec0a463707b1c6
