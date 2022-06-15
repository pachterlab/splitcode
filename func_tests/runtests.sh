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

# Test that program can be run

cmdexec "$splitcode --version"

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


