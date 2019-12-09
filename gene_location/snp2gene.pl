#!/usr/bin/perl -w

$input1=$ARGV[0];
$input2=$ARGV[1];
$output_file=$ARGV[2];

open INPUT1, "$input1";
open INPUT2, "$input2";
open OUTPUT, ">>$output_file";
@raw_data1=<INPUT1>;
@raw_data2=<INPUT2>;

    for ($i=0; $i<=$#raw_data1; $i++)
       {@broken1 = split(/ /, $raw_data1[$i]);
        for ($j=0; $j<=$#raw_data2; $j++)
           {@broken2 = split(/\t/, $raw_data2[$j]);
            if(($broken1[2]>=$broken2[1]) && ($broken1[2]<=$broken2[2])){
               chop($broken2[3]);
               $tmp=$broken1[2]+1;
               print OUTPUT "chrXX $broken1[2] $tmp $broken1[3] $broken1[4] + chrXX.$tmp $broken2[1] $broken2[2] $broken2[3]\n";
               }
           }
       }
close (OUTPUT);
close (INPUT1);

