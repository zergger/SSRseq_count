$| = 1; my $SCRIPT_INFO = "\n"."#"x(30)."
ssr reads count
Version: v1.0 2020-06-21
Contact: 129 ganb
\n"; print $SCRIPT_INFO;
# Perl 系统包
use warnings;
use strict;
use File::Spec;
use Getopt::Long;
use Bio::SeqIO;
use Bio::SearchIO;
use Parallel::ForkManager;
use Data::Dumper; # 导入Data::Dumper模块
use Storable qw(store); # Import 'store' function from Storable

# 常量
use constant SCRIPT_DIR => (File::Spec->splitpath(File::Spec->rel2abs($0)))[1];


# 参数输入
my $ARGV_INFO = join " ", @ARGV;
my ($fastq_dir, $sample_list, $target_fasta, $output_dir, $parallel, $SOFT_FLASH, $SOFT_BLASTN, $SOFT_FASTQ_TO_FASTA, $ploid, $if_help);
GetOptions(
        "fastq_dir|i=s"    => \$fastq_dir,
        "sample_list|s=s"  => \$sample_list,
        "target_seq|t=s"   => \$target_fasta,
        "output_dir|o=s"   => \$output_dir,

        "parallel=s"       => \$parallel,
        "software_flash=s"     => \$SOFT_FLASH,
        "software_blastn=s"    => \$SOFT_BLASTN,
        "software_fastq_to_fasta=s"    => \$SOFT_FASTQ_TO_FASTA,
        "ploid=s"       => \$ploid,
        "help|h"           => \$if_help,
);
die "
Options: [required]

        --fastq_dir/-i                Paired-end fastq file dir
                                      file name format: sample_name_1.fq, sample_name_2.fq
        --target_fasta/-t             pcr targeted fasta seq, seq name should be named as our designed. 
        --output_dir/-o               output dir
        --parallel                    how many samples will be analyzed in parallel, example : 10
        --software_flash              path/to/flash     flash software directory. please download and install by yourself. http://ccb.jhu.edu/software/FLASH/
        --software_blastn             path/to/blastn      blastn software directory. please download and install by yourself.  ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
        --software_fastq_to_fasta     path/to/fastq_to_fasta    fastq_to_fasta software directory. please download and install by yourself. http://hannonlab.cshl.edu/fastx_toolkit/commandline.html

Options: [optional]
        --sample_list/-s          samples that need to be analyzed. If there are multiple samples, separate them with commas. example: sample1,sample2,sample3
                                  we will use all samples in fastq_dir by default.
        --help/-h                 help doc

\n" if (defined $if_help or not defined $fastq_dir or not defined $target_fasta or not defined $output_dir or not defined $parallel or not defined $SOFT_FLASH or not defined $SOFT_BLASTN);

$sample_list = get_sample_list($fastq_dir) if(not defined $sample_list);
# 如果没有设置 ploid 参数，则默认为 2
$ploid = 2 unless defined $ploid;
###################################################################### 初始化
make_dir($output_dir);

my $DATA_TIME = `date +\"\%Y-\%m-\%d \%H:\%M.\%S\"`;
my $RUN_INFO = "
---
Command: perl ".File::Spec->rel2abs($0)." $ARGV_INFO
---
Start time: $DATA_TIME
[SET]  fastq_dir : $fastq_dir
[SET]  sample_list : $sample_list
[SET]  target_fasta : $target_fasta
[SET]  output_dir : $output_dir
[SET]  parallel : $parallel
[SET]  soft_flash : $SOFT_FLASH
[SET]  soft_blastn : $SOFT_BLASTN
[SET]  soft_fastq_to_fasta : $SOFT_FASTQ_TO_FASTA
";
open SAVE, ">>$output_dir/run.info"; 
print SAVE $SCRIPT_INFO.$RUN_INFO; 
close SAVE;
print $RUN_INFO;

die "[Error] lost flash soft\n"          if(is_file_ok($SOFT_FLASH) == 0);
die "[Error] lost blastn soft\n"         if(is_file_ok($SOFT_BLASTN) == 0);
die "[Error] lost fastq_to_fasta soft\n" if(is_file_ok($SOFT_FASTQ_TO_FASTA) == 0);
die "[Error] lost target_fasta\n"        if(is_file_ok($target_fasta) == 0);
die "[Error] lost fastq_dir\n"           if(is_dir_ok($fastq_dir) == 0);

my $blast_plus_dir = `dirname $SOFT_BLASTN`;
   $blast_plus_dir=~s/[\r\n]//g;
my $SOFT_MAKEBLASTDB = "$blast_plus_dir/makeblastdb";
die "[Error] lost makeblastdb soft. please set blastn and makeblastdb at the same directory\n" if(is_file_ok($SOFT_MAKEBLASTDB) == 0);

my @samples = split /,/, $sample_list;

###################################################################### 主程序
check_fastq();  
my $target_motif = check_target_motif();  
my $target_fasta_db = build_blast_index();  

# sample result dir
print "[process] start mapping sample and check fastq str info\n";
my $mapping_dir = "$output_dir/mapping";
make_dir($mapping_dir);

my $pm = Parallel::ForkManager->new($parallel);
foreach my $sample(@samples)
{
    $pm->start() and next;
    sub_run($sample);
    $pm->finish;    
}
$pm->wait_all_children;

print "[process] summary\n";
quality_control();
str_count_summary();


###################################################################### 子程序

sub str_count_summary{
    my %hashSTR;

    # read sample str count
    foreach my $sample(@samples)
    {
        my $str_count_file = "$mapping_dir/$sample/$sample.str_count.txt";
        open STR, $str_count_file;
        my $line1 = <STR>;
           $line1 =~ s/[\r\n]//g;
        my @heads = split /\t/, $line1;

        while(<STR>)
        {
            $_ =~ s/[\r\n]//g;
            next if($_!~/\w/);
            my @datas = split /\t/, $_;
            my %hashTmp = map{ ($heads[$_], $datas[$_])  } (0..$#heads);

            my $target       = $hashTmp{"target"};
            my $aim_str      = $hashTmp{"aim_str"};
            my $str_identify = $hashTmp{"str_identify"};
            my $reads_count  = $hashTmp{"reads_count"};
            # $hashSTR{$target}{$aim_str}{$str_identify}{$sample} = $reads_count;
			
            # Parse aim_str to extract motif
            my ($motif, $start, $end) = split /,/, $aim_str;
            my $motif_length   = length($motif);
            my $str_length     = length($str_identify);
            # my $str_motif_count = $motif_length > 0 ? int($str_length / $motif_length) : 0;
            # # Aggregate counts based on motif and motif count
            # $hashSTR{$target}{$aim_str}{$str_motif_count}{$sample} += $reads_count;
# 如果 motif_length=0，无法计算，直接跳过或做别的处理
if ($motif_length == 0) {
    # 跳过本次循环
    next;
}
    # 计算出 motif count，此时motif count可能不是整数，因为motif的定义不合适，与str实际突变步长不一致。
    my $str_motif_count = $str_length / $motif_length;

    # 先在这种情况下累加 reads_count, **并存储 str_identify**
    $hashSTR{$target}{$aim_str}{$str_motif_count}{$sample}{'Reads'} += $reads_count; # 使用 'Reads' 键来存储 reads_count
    $hashSTR{$target}{$aim_str}{$str_motif_count}{$sample}{'str_identify'} = $str_identify; # **新增：存储 str_identify**
        }
        close STR;
    }

    # # --------------------  **Save %hashSTR to file using Storable**  --------------------
    # my $hashSTR_store_file = "$output_dir/hashSTR.dat"; # Same filename as debug_hashSTR.pl expects
    # store \%hashSTR, $hashSTR_store_file; # Use Storable's store function to serialize and save
    # print "Info: %hashSTR saved to file (Storable format): $hashSTR_store_file\n";
    # # --------------------------------------------------------------------

    # --------------------  动态调整 aim_str (Motif)  --------------------
    # **检测 $hashSTR 中非整数 $str_motif_count 情况**

    my %total_reads_aim_str; # 用于存储每个 aim_str 在 *每个样本中* 的总 reads 数量 (修改后的定义)

    # **首先计算每个 aim_str 在 *每个样本中* 的总 reads 数量 (修改后的计算逻辑)**
    foreach my $target_key (keys %hashSTR) {
        foreach my $aim_str_key (keys %{$hashSTR{$target_key}}) {
            foreach my $sample (@samples) { # **遍历 @samples 列表**
                my $current_aim_str_total_reads = 0; # 初始化当前 aim_str 在 *当前样本* 中的总 reads 计数器
                foreach my $str_motif_count_key (keys %{$hashSTR{$target_key}{$aim_str_key}}) {
                    if (exists $hashSTR{$target_key}{$aim_str_key}{$str_motif_count_key}{$sample}{'Reads'}) {
                        $current_aim_str_total_reads += $hashSTR{$target_key}{$aim_str_key}{$str_motif_count_key}{$sample}{'Reads'}; # 累加 reads count
                    }
                }
                $total_reads_aim_str{$aim_str_key}{$sample} = $current_aim_str_total_reads; # **存储 aim_str 在 *当前样本* 中的总 reads 数量, 使用嵌套 hash 结构**
            }
        }
    }

    my %hash_non_integer_motif_count; # 用于记录出现非整数motif count的 aim_str
    foreach my $target_key (keys %hashSTR) {
        foreach my $aim_str_key (keys %{$hashSTR{$target_key}}) {
            foreach my $str_motif_count_key (keys %{$hashSTR{$target_key}{$aim_str_key}}) {
                foreach my $sample (@samples) { # **遍历 @samples 列表，确保针对每个样本进行判断**
                    if (exists $hashSTR{$target_key}{$aim_str_key}{$str_motif_count_key}{$sample}{'Reads'}) {
                        my $reads_count = $hashSTR{$target_key}{$aim_str_key}{$str_motif_count_key}{$sample}{'Reads'};
                        my $aim_str_total_reads = $total_reads_aim_str{$aim_str_key}{$sample}; # **获取当前 aim_str 在 *当前样本* 中的总 reads 数量**
# print $aim_str_key."====".$str_motif_count_key."====".$reads_count / $aim_str_total_reads."\n";
                        # **判断条件：reads_count 是否超过总 reads 数量的 (1/$ploid-0.1)**
                        if ($aim_str_total_reads > 0 && ($reads_count / $aim_str_total_reads) >= (1/$ploid-0.1)) {
                            # 如果满足条件，则执行非整数 motif_count 检测逻辑
                            if ($str_motif_count_key ne int($str_motif_count_key)) { # 非整数判断
                                # 发现非整数 motif count
                                # print "Warning: Non-integer motif_count found for target: $target_key, aim_str: $aim_str_key, motif_count: $str_motif_count_key\n";
                                $hash_non_integer_motif_count{$aim_str_key} = $str_motif_count_key; # 存储 str_motif_count_key 本身作为值
                            }
                        }
                    }
                }
            }
        }
    }

    # print "Warning: Found aim_str with non-integer motif counts. Will attempt to adjust aim_str.\n";
    my %hash_aim_str_to_new_motif; # 用于存储旧 aim_str 到新 motif 的映射

    # **修改 aim_str (Motif) 定义**
    foreach my $aim_str_key (keys %hash_non_integer_motif_count) {
        # 从 aim_str 中解析出 motif, start, end
        my ($motif, $start, $end) = split /,/, $aim_str_key;
        my $motif_length = length($motif);

        # **直接从 %hash_non_integer_motif_count 哈希中获取 str_motif_count_key**
        my $non_integer_motif_count_key = $hash_non_integer_motif_count{$aim_str_key};

        #  新的 motif 取原 motif 的前 motif_length/2 或 motif_length/3 长度的碱基
        if (defined $non_integer_motif_count_key && $motif_length > 1) { # 确保成功获取到非整数 motif_count 且 motif 长度大于 1
            my $decimal_part = $non_integer_motif_count_key - int($non_integer_motif_count_key); # 计算小数部分

            my $new_motif_length;
            if ($motif_length == 3 || $motif_length == 6) { # **针对 motif_length 为 3 或 6 的情况进行特殊处理**
                if ((abs($decimal_part - 0.66666) < 0.05) || (abs($decimal_part - 0.33333) < 0.05)) {
                    if ($motif_length == 3) {
                        $new_motif_length = $motif_length; #  motif_length 为 3 时，不改变 motif 长度
                    } elsif ($motif_length == 6) {
                        $new_motif_length = int($motif_length / 3); # motif_length 为 6 时，motif 长度调整为 1/3
                    }
                } else {
                    $new_motif_length = $motif_length; # 小数部分不符合 0.666 或 0.333， 默认不改变 motif 长度
                }
            } elsif (abs($decimal_part - 0.5) < 0.05) {
                $new_motif_length = int($motif_length / 2); # 减半
            } else {
                $new_motif_length = $motif_length; # **对于其他 motif_length，不改变 motif 长度**
            }

            if ($new_motif_length > 0) {
                my $new_motif = substr($motif, 0, $new_motif_length);
                my $new_aim_str = join(",", $new_motif, $start, $end);
                $hash_aim_str_to_new_motif{$aim_str_key} = $new_aim_str;
            }

        } else {
            # print "Warning: aim_str '$aim_str_key' motif length is too short ($motif_length) to adjust. Skipping.\n";
            #  如果 motif 长度太短，可以跳过
        }
    }

    my %hashSTR_new; #  新的 hashSTR

    # **调整 $hashSTR 基于新的 aim_str 定义**
    foreach my $target_key (keys %hashSTR) {
        foreach my $aim_str_key (keys %{$hashSTR{$target_key}}) {
            foreach my $str_motif_count_key (keys %{$hashSTR{$target_key}{$aim_str_key}}) {
                foreach my $sample (@samples) {  #  注意这里需要遍历 @samples 获取样本名
                    if (exists $hashSTR{$target_key}{$aim_str_key}{$str_motif_count_key}{$sample}) { # 检查样本是否存在数据
                        my $reads_count = $hashSTR{$target_key}{$aim_str_key}{$str_motif_count_key}{$sample}{'Reads'};
                        my $str_identify = $hashSTR{$target_key}{$aim_str_key}{$str_motif_count_key}{$sample}{'str_identify'};

                        my $current_aim_str_key = $aim_str_key; #  当前处理的 aim_str
                        my $current_motif;
                        my $current_motif_length;

                        if (exists $hash_aim_str_to_new_motif{$aim_str_key}) {
                            # 如果 aim_str 需要调整，则使用新的 aim_str 和 motif
                            $current_aim_str_key = $hash_aim_str_to_new_motif{$aim_str_key};
                            ($current_motif, my $start, my $end) = split /,/, $current_aim_str_key;
                            $current_motif_length = length($current_motif);
                        } else {
                            # aim_str 不需要调整，使用原始的 motif
                            ($current_motif, my $start, my $end) = split /,/, $current_aim_str_key;
                            $current_motif_length = length($current_motif);
                        }

                        #  重新计算 motif_length 和 str_motif_count (使用新的或原始的 motif 定义)
                        my $current_str_length = length($str_identify);
                        my $current_str_motif_count; # 先声明 current_str_motif_count 变量

                        # 判断 $current_str_length 是否能被 $current_motif_length 整除，并计算 $current_str_motif_count
                        if (exists $hash_aim_str_to_new_motif{$aim_str_key}) {
                            # 如果 $hash_aim_str_to_new_motif{$aim_str_key} 存在，则直接取整计算 motif count
                            if ($current_motif_length > 0 && defined $reads_count) { # 仅需判断 $current_motif_length > 0 和 $reads_count 是否定义
                                $current_str_motif_count = int($current_str_length / $current_motif_length); # 直接取整
                                #  累加 reads_count 到新的 hashSTR 结构中，并存储 str_identify
                                $hashSTR_new{$target_key}{$current_aim_str_key}{$current_str_motif_count}{$sample}{'Reads'} += $reads_count;
                                $hashSTR_new{$target_key}{$current_aim_str_key}{$current_str_motif_count}{$sample}{'str_identify'} = $str_identify;
                            }
                        } else {
                            # 否则，仍然使用之前的整除判断逻辑
                            if ($current_motif_length > 0 && defined $reads_count && $current_str_length % $current_motif_length == 0) { # 添加条件判断
                                # 如果能整除，则计算 $current_str_motif_count
                                $current_str_motif_count = $current_str_length / $current_motif_length; #  重新计算 motif count (现在只在整除时计算)
                                #  累加 reads_count 到新的 hashSTR 结构中
                                $hashSTR_new{$target_key}{$current_aim_str_key}{$current_str_motif_count}{$sample}{'Reads'} += $reads_count;
                                $hashSTR_new{$target_key}{$current_aim_str_key}{$current_str_motif_count}{$sample}{'str_identify'} = $str_identify;
                            } else {
                                # 如果不能整除，则跳过当前数据条目
                                next;
                            }
                        }

                    }
                }
            }
        }
    }

    #  用新的 hashSTR 替换旧的 hashSTR
    %hashSTR = %hashSTR_new;
    # print "Info: hashSTR adjusted based on new aim_str definitions.\n";

    # output
    my $str_count_summary = "$output_dir/str_count.txt";
    open OUT, ">$str_count_summary";

    print OUT "Target\tAimSTR\tSTR\t" . (join "\t", @samples) . "\n";
    foreach my $target(sort keys %hashSTR)
    {
        foreach my $aim_str(sort keys %{$hashSTR{$target}})
        {   
            my ($motif, $start, $end) = split /,/, $aim_str;
            my $motif_length = length($motif);
            # foreach my $str_identify (sort keys %{$hashSTR{$target}{$aim_str}})
			foreach my $str_motif_count (sort { $a <=> $b } keys %{ $hashSTR{$target}{$aim_str} })
            {
                # my $str_motif_count = int(length($str_identify) / $motif_length);
                my @values = ($target, $aim_str, "$motif($str_motif_count)");
                foreach my $sample(@samples)
                {
                    # my $reads_count = exists $hashSTR{$target}{$aim_str}{$str_identify}{$sample} ? $hashSTR{$target}{$aim_str}{$str_identify}{$sample} : "";
					my $reads_count = exists $hashSTR{$target}{$aim_str}{$str_motif_count}{$sample} 
                                      ? $hashSTR{$target}{$aim_str}{$str_motif_count}{$sample}{'Reads'} # 修改后的代码，访问 'Reads' 键
                                      : "";
                    push @values, $reads_count;
                }
                print OUT (join "\t", @values) . "\n";
            }
        }
    }
    close OUT;

    print "STR Count: $str_count_summary\n";
    
}

sub quality_control{
    my %hashQC;
    my %hashQCT;

    foreach my $sample(@samples)
    {
        my $flash_log         = "$mapping_dir/$sample/flash.log.txt";
        my $statistics_target = "$mapping_dir/$sample/$sample.filter_ok.txt";

        open FLASH, $flash_log;
        while(<FLASH>)
        {
            if($_=~/Total pairs:\s+(\d+)/)
            {
                $hashQC{$sample}{'RawReads'} = $1;
                next;
            }
            if($_=~/Combined pairs:\s+(\d+)/)
            {
                $hashQC{$sample}{'MergeReads'} = $1;
                next;
            }
        }
        close FLASH;

        open TARGET, $statistics_target;
        while(<TARGET>)
        {
            $_=~s/[\r\n]//g;
            my ($target, $count) = split /\t/, $_;
            $hashQCT{$target}{$sample} = $count;
            $hashQC{$sample}{'FilterReads'} += $count;
        }
        close TARGET;

        $hashQC{$sample}{'Sample'} = $sample;
        $hashQC{$sample}{'MergePerc'} = sprintf "%0.2f", $hashQC{$sample}{'MergeReads'} / $hashQC{$sample}{'RawReads'};
        $hashQC{$sample}{'FilterPerc'} = sprintf "%0.2f", $hashQC{$sample}{'FilterReads'} / $hashQC{$sample}{'MergeReads'};

    }

    # output result
    my $qc = "$output_dir/qc.txt";
    my @heads_qc = qw(Sample RawReads MergeReads MergePerc FilterReads FilterPerc);
    open QC, ">$qc";
    print QC (join "\t", @heads_qc) . "\n";
    foreach my $sample(@samples)
    {
        my @values = map{ $hashQC{$sample}{$_} } @heads_qc;
        print QC (join "\t", @values) . "\n";
    }

    my $qct = "$output_dir/qc_target.txt";
    my @targets = sort keys %hashQCT;

    open QCT, ">$qct";
    print QCT "Sample\t" . (join "\t", @targets) . "\n";
    foreach my $sample(@samples)
    {
        my @counts = map{ $hashQCT{$_}{$sample} } @targets;
        print QCT "$sample\t" . (join "\t", @counts) . "\n";
    }
    close QCT;

    print "QC: $qc\n";
    print "QC Target: $qct\n";
}

 

sub sub_run{
    my $sample   = shift @_;
    my $fastq_r1 = "$fastq_dir/$sample\_1.fq.gz";
    my $fastq_r2 = "$fastq_dir/$sample\_2.fq.gz";

    my $sample_dir = "$mapping_dir/$sample";
    make_dir($sample_dir);

    # Merge
    print "1.Merge Fastq : $sample\n";
    # system "$SOFT_FLASH -t 4 --allow-outies -m 15 -M 150 -x 0.1 -d $sample_dir -o $sample $fastq_r1 $fastq_r2 > $sample_dir/flash.log.txt 2>&1"; 
	# M提高，能一定程度上提高短片段的合并成功数量，提高体系富集效率大约1%，但对总体来说几乎无影响
	system "$SOFT_FLASH -t 4 --allow-outies -M 150 -d $sample_dir -o $sample $fastq_r1 $fastq_r2 > $sample_dir/flash.log.txt 2>&1";

    # filter
    print "2.Filter Fastq : $sample \n";
    my $fastq_merge_mapped    = "$sample_dir/$sample.merge.mapped.fastq";
    my $merge_reads_belong    = "$sample_dir/$sample.merge.reads.belong";
    my %hashTargetFasta  = read_fasta($target_fasta);
	# print Dumper(\%hashTargetFasta);
	# exit;

	my $not_combined_1 = "$sample_dir/$sample.notCombined_1.fastq";
	my $not_combined_2 = "$sample_dir/$sample.notCombined_2.fastq";
	my $no_overlap_merged = "$sample_dir/$sample.notCombined_merged.fastq";
	# 调用函数进行无重叠合并
	merge_no_overlap_files($not_combined_1, $not_combined_2, $no_overlap_merged);

    merge_fastq_analysis($target_fasta_db, \%hashTargetFasta, $sample_dir, $sample); # 合并成功reads过滤

    my $final_fastq = "$sample_dir/$sample.fastq";
    system("cat $fastq_merge_mapped > $final_fastq");

    # fastq比对
    print "3.blast final Fastq : $sample \n";
    # my $final_blast = "$sample_dir/$sample.blast.gz";
	my $final_blast = "$sample_dir/$sample.blast";
    blast_final($final_fastq, $merge_reads_belong, $target_fasta_db, $final_blast, $sample_dir);

    # str reads count 
    print "4.get str count : $sample \n";
    str_count($merge_reads_belong, $final_blast, $sample_dir, $sample);

}


sub str_count{
    my $reads_belong = shift @_;
    my $blast_out    = shift @_;
    my $sample_dir   = shift @_;
    my $sample       = shift @_;


    my %hashMotif       = read_motif($target_motif);
    my %hashReadsBelong = read_read_belong($reads_belong, ""); # 经过上一步的质控，提取reads属于哪一个参考片段
    
    ### 参照序列STR信息，只需要起始与终止位置
    my %hashResult;   # 记录str信息
    my %hashBlastSeq; # 记录比对序列，用于后期核查
    my %hashFinishReads; # reads 已经处理完毕
    
    # open my $fh, "gzip -cd $blast_out |";
	open my $fh, '<', $blast_out;
    my $blast_in = Bio::SearchIO->new(-fh => $fh, -format => 'blast');
    while( my $r = $blast_in->next_result ) {
        while( my $h = $r->next_hit ) {  
            while( my $hsp = $h->next_hsp ) {
				# print Dumper($h);
				# print "=======我是分割线=======\n";
				# print Dumper($hsp);
				# exit;
                my $hit_name   = $h->name;       # 参照序列名称
                my $query_name = $r->query_name ;# reads序列名称
                   $hit_name   =~ s/\s+//g;
                   $query_name =~ s/\s+//g;
				   $hit_name   =~ s/\:/\|/g;
				# print "$hit_name===$query_name\n";
				# exit;
                last if(exists $hashFinishReads{$query_name} ); # 该reads已经分析完毕
                last if($hit_name ne $hashReadsBelong{$query_name}); # 该reads质控时对应的参考片段不是当前的参考片段

                my $hit_string   = uc($hsp->hit_string);  # 比对序列,参考序列
                my $qerry_string = uc($hsp->query_string); # 比对序列，reads
                my ($hit_start, $hit_end) = $hsp->range('hit') ; # 参考序列比对范围，永远是从小到大的顺序，即使是Minus
                my $fangxiang = $hsp->strand('hit') ; # 获取参照序列的方向 
                # 如果参照序列是反向匹配的，则需要反转,必须保证参考序列是正向的
                if($fangxiang == -1){ 
                   $qerry_string = reverse $qerry_string;
                   $qerry_string =~ tr/ATCGatcg/TAGCtagc/;
                   $hit_string   = reverse $hit_string; 
                   $hit_string   =~ tr/ATCGatcg/TAGCtagc/;                                        
                }

                # 判定该reads的STR序列
                foreach my $aim_str(keys %{$hashMotif{$hit_name}})# 对需要在该片段上判定的目标STR循环识别
                {   
                    # 获取参考序列上STR情况
                    my $motif         = $hashMotif{$hit_name}{$aim_str}{'motif'};
                    my $motif_length  = $hashMotif{$hit_name}{$aim_str}{'motif_length'};
                    my $str_seq_start = $hashMotif{$hit_name}{$aim_str}{'str_seq_start'};
                    my $str_seq_end   = $hashMotif{$hit_name}{$aim_str}{'str_seq_end'};
                    
                    # 比对区域没有包含目标STR区域
                    next if($str_seq_start < $hit_start or $str_seq_end > $hit_end);
                    
                    # # 获取参照STR两端各扩展3bp，因为当indel出现时，可能会造成原始序列STR变长
                    # my $extend = 3;
                    # $str_seq_start = $str_seq_start - $extend;
                    # $str_seq_end   = $str_seq_end + $extend;

                    # (1) 第一步，循环遍历hit_string, 根据参考位置，提取目标STR的序列
                    my $get_seq     = ""; # 实际reads在目标STR区域的序列
                    my $before_seq  = ""; # 目标区域前面的序列
                    my $after_seq   = ""; # 目标区域后面的序列
                    my $target_pos  = $hit_start; # 参考序列上的绝对位置

# 统计缺口数量及其位置
my @gap_positions = ();
for (my $i = 0; $i < length($hit_string); $i++) {
    my $hit_base = substr($hit_string, $i, 1);
    if ($hit_base eq '-') {
        push @gap_positions, $i + 1; # 修改为1-based位置
    }
}
# 处理缺口并调整目标区域
if (@gap_positions) {
    foreach my $gap_pos (@gap_positions) {
        # 将gap_pos转换为参考序列的位置
        my $gap_ref_pos = $hit_start + $gap_pos;
		if ($gap_ref_pos > $str_seq_start) {
            # 缺口在目标区域之后，增加$str_seq_end
            # 每个缺口增加1
            $str_seq_end += 1;
        }
    }
    # print "目标区域已根据缺口位置调整为：$str_seq_start 到 $str_seq_end\n";
} else {
    # print "未检测到缺口，目标区域保持不变：$str_seq_start 到 $str_seq_end\n";
}
# 确保$str_seq_start和$str_seq_end的有效性
$str_seq_start = 1 if $str_seq_start < 1;
$str_seq_end = length($hit_string) if $str_seq_end > length($hit_string);

                    for(my $hit_pos = 0; $hit_pos < length($hit_string); $hit_pos++)
                    {
                        my $hit_base   = substr($hit_string, $hit_pos, 1);
                        my $query_base = substr($qerry_string, $hit_pos, 1);
                        $get_seq    = $get_seq.$query_base    if($target_pos >= $str_seq_start and $target_pos <= $str_seq_end);
                        $before_seq = $before_seq.$query_base if($target_pos < $str_seq_start);
                        $after_seq  = $after_seq.$query_base  if($target_pos > $str_seq_end);
                        # $target_pos++ if($hit_base=~/[ATCG]/i);# 只有当hit_base是碱基的时候，绝对位置向后移动一位，因为存在插入缺失，即：字符“-”
						$target_pos++;
                    }
					# print Dumper($query_name);
					# print Dumper($hit_name);
					# print Dumper($hit_string);
					# print Dumper($get_seq);
					# print Dumper($before_seq);
					# print Dumper($after_seq);
					# exit;
                    # (2) 第二步，判定提取的序列上STR信息
                    # 去掉插入缺失字符
                    $get_seq    =~ s/[-]//g;
                    $before_seq =~ s/[-]//g;
                    $after_seq  =~ s/[-]//g;
                    my ($str_seq, $start, $end, $str_seq_length) = get_seq_str($query_name, $get_seq, $motif, 1, 'Force'); # 获取该片段最大的STR序列,最少重复1次

                    # # (3) 第三步，从目标STR区域向两侧扩展，每次左右各扩展time * motifLength 长度，因为有可能前面的extend扩展的不足
                    # my $time     = 0; # 左右扩展次数，
                    # my $continue = 1; # 控制是否继续循环
                    # while($continue){
                       # $time++;
                       # my $cut_length    = $motif_length * $time; # 本次前后扩展长度
                       # my $before_cutseq = (length($before_seq) <= $cut_length) ? $before_seq : substr($before_seq, length($before_seq) - $cut_length, $cut_length); # getSeq前面扩展的序列
                       # my $after_cutseq  = (length($after_seq) <= $cut_length)  ? $after_seq : substr($after_seq, 0, $cut_length); # getSeq后面扩展的序列
                       # my $combine_seq   = $before_cutseq.$get_seq.$after_cutseq;                  
                       # my ($str_seq_tmp, $start_tmp, $end_tmp, $str_seq_length_tmp) = get_seq_str($query_name, $combine_seq, $motif, 1, 'Force'); # 获取扩展后该片段最大的STR序列

                       # # 扩展后变长
                       # if($str_seq_length < $str_seq_length_tmp){
                          # $str_seq_length = $str_seq_length_tmp;
                          # $str_seq        = $str_seq_tmp;
                       # }else{ # 扩展并不会使STR变长，离开循环
                          # $continue = 0;
                       # }
                    # }
                    # (4) 第四步，结果保存
                    if($str_seq =~ /[agtc]/ig){
                        $hashFinishReads{$query_name}++; # 记住当前reads名称，表明其已分析完毕，避免同源导致多次处理
                        $hashResult{$hit_name}{$aim_str}{$str_seq}++; # 参考片段上记录一个str
                        $hashBlastSeq{$hit_name}{$aim_str}{$str_seq}{$query_name}{'querry'} = $qerry_string; # 详细比对序列，用于后续输出，检查
                        $hashBlastSeq{$hit_name}{$aim_str}{$str_seq}{$query_name}{'hit'}    = $hit_string;
                    }
                }                   
            }
        }
    }
    close $fh;


    # 结果输出
    open STRCOUNT, ">$sample_dir/$sample.str_count.txt"; # 每种STR的reads数量
    open STRCOUNT_DETAIL, "| gzip > $sample_dir/$sample.str_count.detail.xls.gz"; # 每条reads分型的STR详情，用于检查是否分型有问题。压缩输出，减少空间浪费
    open READ_LENGTH, ">$sample_dir/$sample.reads.length.count.txt"; # 每个片段上的reads每种长度reads的数量统计

    print STRCOUNT "target\taim_str\tstr_referrence\tstr_identify\treads_count\n";
    print STRCOUNT_DETAIL "target\taim_str\tstr_identify\treads_count\treads_name\tseq\n";
    print READ_LENGTH "target\taim_str\tlength\tcount\n";
   
    foreach my $target(sort keys %hashResult)
    {
        foreach my $aim_str(sort keys %{$hashResult{$target}})
        {   
            my %hashLength;
            foreach my $str_identify(sort keys %{$hashResult{$target}{$aim_str}})
            {
                print STRCOUNT "$target\t$aim_str\t$hashMotif{$target}{$aim_str}{'str_seq'}\t$str_identify\t$hashResult{$target}{$aim_str}{$str_identify}\n";
                foreach my $querryname(sort keys %{$hashBlastSeq{$target}{$aim_str}{$str_identify}}){
                    print STRCOUNT_DETAIL "$target\t$aim_str\t$str_identify\t$hashResult{$target}{$aim_str}{$str_identify}\t$querryname\_QUERRY\t$hashBlastSeq{$target}{$aim_str}{$str_identify}{$querryname}{'querry'}\n";
                    print STRCOUNT_DETAIL "$target\t$aim_str\t$str_identify\t$hashResult{$target}{$aim_str}{$str_identify}\t$querryname\_HIT\t$hashBlastSeq{$target}{$aim_str}{$str_identify}{$querryname}{'hit'}\n";            
                    
                    # 记录序列长度
                    my $seq = $hashBlastSeq{$target}{$aim_str}{$str_identify}{$querryname}{'querry'};
                       $seq =~ s/-//g;
                    my $length = length($seq);
                    $hashLength{$length}++;
                }
            }
            # 输出长度信息
            my $motif_length = $hashMotif{$target}{$aim_str}{'motif_length'};
            my @lengths      = sort {$a <=> $b} keys %hashLength;
            my $min_length   = $lengths[0] - $motif_length;
            my $max_length   = $lengths[$#lengths] + $motif_length;
            $hashLength{$min_length}++;
            $hashLength{$max_length}++;# 最大、最小长度上下扩展1个motif，使图形更美观
            map{ print READ_LENGTH "$target\t$aim_str\t$_\t$hashLength{$_}\n"; } sort {$a <=> $b} keys %hashLength; # 输出reads长度的数量            
        }
    }
    close STRCOUNT;
    close STRCOUNT_DETAIL;
    close READ_LENGTH;
}



# 最终比对，按片段拆分执行，以节省CPU/存储消耗
sub blast_final{
    my $final_fastq        = shift @_; # 要处理的fastq
    my $reads_belong_merge = shift @_; # 每条reads所属target
    my $target_fasta_db    = shift @_; # target数据库
    my $final_blast        = shift @_; # 最终输出文件
    my $output_dir         = shift @_; 

    # 临时目录，用后即删
    my $tmp_dir = "$output_dir/fastq_split_blast"; 
    make_dir($tmp_dir);


    my %hashTargetFasta = read_fasta($target_fasta_db); # 参考序列名称获取 {'target_name'}= seq
    my %hashReadsBelong = read_read_belong($reads_belong_merge, ""); # reads所属片段信息获取 {'reads_name'} = target_name

    # 创建每个片段fa句柄
    my %hashHandle; 
    foreach my $target_name(keys %hashTargetFasta)
    {
        my ($target_name_prefix) = split /[|]/, $target_name;  # 取片段前缀作为文件名
        my $tmp_fa    = "$tmp_dir/$target_name_prefix.fa";         # 属于当前片段的reads序列
        my $tmp_limit = "$tmp_dir/$target_name_prefix.limit.txt";  # 当前片段的完整名称，用于blastn的输入
		# 替换 $target_name 中的 | 为 :
		my $target_name_limit = $target_name;
		$target_name_limit =~ s/\|/:/g;
        # system("echo '$target_name\n' > $tmp_limit");
        system("echo '$target_name_limit\n' > $tmp_limit");
        open $hashHandle{$target_name}, ">$tmp_fa";
    }

    # final_fastq 拆分到每个target片段句柄里
    open FASTQ, $final_fastq;
    while(my $line1 = <FASTQ>)
    {
        my $line2 = <FASTQ>;
        my $line3 = <FASTQ>;
        my $line4 = <FASTQ>;

           $line1=~ s/[\r\n]//g;
           $line1=(split /\s/,$line1)[0];
           $line1=~ s/^\@//g;
        my $handle = $hashHandle{$hashReadsBelong{$line1}};
        print $handle ">$line1\n$line2";
    }
    close FASTQ;  

    # 关闭句柄
    foreach my $target_name(keys %hashHandle)
    {
        close $hashHandle{$target_name};
    }
    %hashReadsBelong = (); # 清空，释放内存

    # 比对
    my @blast_results; # 所有比对结果文件
    foreach my $target_name(keys %hashTargetFasta)
    {
        my ($target_name_prefix) = split /[\|]/, $target_name;
        my $tmp_fa    = "$tmp_dir/$target_name_prefix.fa";         # 属于当前片段的reads序列
        my $tmp_limit = "$tmp_dir/$target_name_prefix.limit.txt";  # 当前片段的名称，用于blastn的输入
		my $tmp_limit_in = "$tmp_dir/$target_name_prefix.limit.in.txt";  # 当前片段的名称，用于blastn的输入
        # my $tmp_blast = "$tmp_dir/$target_name_prefix.blast.gz";
        my $tmp_blast = "$tmp_dir/$target_name_prefix.blast";
        next if(is_file_ok($tmp_fa) == 0); # 没有数据，不必再做
        # system("$SOFT_BLASTN -task blastn -query $tmp_fa -evalue 0.00001 -db $target_fasta_db -num_threads 4 -dust no -seqidlist $tmp_limit | gzip > $tmp_blast"); # 压缩，减少存储消耗
		system("blastdb_aliastool -seqid_file_in $tmp_limit -seqid_file_out $tmp_limit_in"); # 压缩，减少存储消耗
        system("$SOFT_BLASTN -task blastn -query $tmp_fa -evalue 0.00001 -db $target_fasta_db -num_threads 4 -dust no -seqidlist $tmp_limit_in > $tmp_blast"); # 压缩，减少存储消耗
        push @blast_results, $tmp_blast;
    }   

    # 合并
    system("cat @blast_results > $final_blast");
    # 清空中间文件
    system("rm -r $tmp_dir");
}

# 提取reads所属信息
sub read_read_belong{
    my $reads_belong = shift @_;
    my $type         = shift @_;

    my %hashReadsBelong;
    open BELONG, $reads_belong;
    while(<BELONG>)
    {
        $_=~s/[\r\n]//g;
        my ($reads_name, $target_name) = split /\t/, $_;        
        if($type eq "FASTQ_TYPE") # 原始fastq命名格式
        {
            $reads_name=~s/:MERGE$//;
            $reads_name=~s/:R\d$//;
            $reads_name="\@$reads_name";
        }
        elsif($type eq 'READS_NAME_TYPE')
        {
            $reads_name=~s/:MERGE$//;
            $reads_name=~s/:R\d$//;            
        }
        $hashReadsBelong{$reads_name} = $target_name;
    }
    close BELONG;
    return %hashReadsBelong;
}


# 合并成功reads过滤
sub merge_fastq_analysis{
    my $target_fasta     = shift @_;
    my $hashTargetFasta  = shift @_;
    my $sample_dir       = shift @_;
    my $sample           = shift @_;

    # notCombined 也加进来
    my $merge_fastq = "$sample_dir/$sample.extendedFrags.fastq";
	my $no_overlap_merged = "$sample_dir/$sample.notCombined_merged.fastq";
	system("cat $no_overlap_merged >> $merge_fastq");
    # system("cat $sample_dir/$sample.notCombined_1.fastq >> $merge_fastq");
    # system("cat $sample_dir/$sample.notCombined_2.fastq >> $merge_fastq");

    # fastq 转 fasta
    # my $merge_fastq = "$sample_dir/$sample.extendedFrags.fastq";
    my $merge_fasta = "$sample_dir/$sample.extendedFrags.fa";   
    system "$SOFT_FASTQ_TO_FASTA -Q 33 -n -i $merge_fastq -o $merge_fasta";

    # blast
    my $merge_blast = "$sample_dir/$sample.extendedFrags.blast";
    system("$SOFT_BLASTN -task blastn -query $merge_fasta -outfmt 6 -evalue 0.00001 -db $target_fasta -out $merge_blast -num_threads 4 -max_target_seqs 5 -dust no");
    
    # 筛选
    my %hashmergeFastq = read_fastq($merge_fastq);
    my %hashmergeBlast = read_blast(\%hashmergeFastq, $hashTargetFasta, $merge_blast, 'Merge');
    
    my $fastq_mapped    = "$sample_dir/$sample.merge.mapped.fastq";
    my $fastq_un_mapped = "$sample_dir/$sample.merge.unmapped.fastq";
    my $reads_belong    = "$sample_dir/$sample.merge.reads.belong";
    open MAPPED, ">$fastq_mapped";
    open UNMAPPED, ">$fastq_un_mapped";
    open BELONG, ">$reads_belong";
 
    my %hashFilterOK; # qualified reads count in each target
    foreach my $reads_name(keys %hashmergeFastq){
        if(exists $hashmergeBlast{$reads_name}){
            print MAPPED "\@$reads_name:MERGE\n";
            print MAPPED "$hashmergeFastq{$reads_name}{2}\n";
            print MAPPED "$hashmergeFastq{$reads_name}{3}\n";
            print MAPPED "$hashmergeFastq{$reads_name}{4}\n";
            print BELONG "$reads_name:MERGE\t$hashmergeBlast{$reads_name}\n";

            $hashFilterOK{$hashmergeBlast{$reads_name}}++;

        }else{
            print UNMAPPED "\@$reads_name\n";
            print UNMAPPED "$hashmergeFastq{$reads_name}{2}\n";
            print UNMAPPED "$hashmergeFastq{$reads_name}{3}\n";
            print UNMAPPED "$hashmergeFastq{$reads_name}{4}\n";
        }
    }
    close MAPPED;
    close UNMAPPED;
    close BELONG;

    my $statistics_target = "$sample_dir/$sample.filter_ok.txt";
    open STAT, ">$statistics_target";
    foreach my $target_name(sort keys %$hashTargetFasta)
    {
        my $filter_ok = exists $hashFilterOK{$target_name} ? $hashFilterOK{$target_name} : 0;
        print STAT "$target_name\t$filter_ok\n";
    }
    close STAT;
}

# 读取fastq文件
sub read_fastq{
    my $fastq = shift @_;
    my %hashFastq;
    open FASTQ,$fastq;
    while(my $line1 = <FASTQ>){
        my $line2 = <FASTQ>;
        my $line3 = <FASTQ>;
        my $line4 = <FASTQ>;
        $line1=~ s/[\r\n]//g;
        $line2=~ s/[\r\n]//g;
        $line3=~ s/[\r\n]//g;
        $line4=~ s/[\r\n]//g;
        my $name=(split /\s/,$line1)[0];
           $name=~ s/^\@//;
        $hashFastq{$name}{1} = $line1;
        $hashFastq{$name}{2} = $line2;
        $hashFastq{$name}{3} = $line3;
        $hashFastq{$name}{4} = $line4;
    }
    close FASTQ;
    return %hashFastq;
}

# 读取blast结果
sub read_blast{
    my $hashFastq       = shift @_; 
    my $hashTargetFasta = shift @_;
    my $blast           = shift @_;
    my $type            = shift @_;###single表示对merge的reads进行比对。pair表示对没有merge上的R1或R2的reads进行比对。
    
    # 读取blast比对信息，并进行基础过滤
    my %hashBlast;
    open BLAST,$blast;
    while(<BLAST>){
		chomp; # 去除行末的换行符
		# 替换 : 为 |
		s/:/\|/g;
        # $_=~ s/[\r\n]//g;
        my ($reads_name, $target_name, $identify_perc, $identify_len, $mismatch, $gap, $reads_start, $reads_end, $target_start, $target_end, $evalue, $score) = split /\t/, $_;
        my @array=split(/\t/,$_);
        my $reads_length  = length( $hashFastq->{$reads_name}{2} );
        my $target_length = length( $hashTargetFasta->{$target_name} );
		# print Dumper($reads_name);
		# print Dumper($target_start);
		# print Dumper($target_end);
		# exit;
        ($target_start, $target_end) = ($target_end, $target_start) if($target_end < $target_start);
 
        my $judge_value = $identify_len; # 默认用匹配长度作为后续候选的筛选规则
        #
        # 情况一：Reads覆盖整个目标区域
        # 条件，覆盖超过80%的区域
        # 条件，覆盖起始位置小于开头加5，结束位置大于末尾减5 #这个没必要
        #
        # if($type eq 'Merge' and $identify_len > $target_length * 0.8 and $target_start <= 5 and $target_end > $target_length - 5){
		if($type eq 'Merge' and $identify_len > $target_length * 0.8){
            $hashBlast{$reads_name}{$judge_value} = $target_name;
        }
    }
    close BLAST;

    # 进一步筛选
    my %hashBest=();

    foreach my $title(keys %hashBlast){
        foreach my $len (sort {$b <=> $a} keys %{$hashBlast{$title}}){
            $hashBest{$title} = $hashBlast{$title}{$len};
            last;
        }
    }       

    return %hashBest;
}



sub check_fastq{
    print "[process] Check fastq \n";
    foreach my $sample(@samples)
    {
        my $fastq_r1 = "$fastq_dir/$sample\_1.fq.gz";
        my $fastq_r2 = "$fastq_dir/$sample\_2.fq.gz";
        die "[Error] lost fastq: $fastq_r1 $fastq_r2\n" if(is_file_ok($fastq_r1, $fastq_r2) == 0);
    }
}

sub check_target_motif{
    my $hashConfig = shift @_;

    print "[process] Check target motif \n";
    my %hashTarget   = read_fasta($target_fasta);

    # （1）提取motif详细信息
    my %hashResult; # 记录该片段需要查找的所有motif，并判断是否重复
    my $count = 0;
    foreach my $target_full_name(sort keys %hashTarget)
    {
        my $target_seq        = $hashTarget{$target_full_name};
        my $target_seq_length = length($hashTarget{$target_full_name});# 参考序列长度
        my ($target_name, $motif_info_list) = split /\|/, $target_full_name, 2;
        
        foreach my $motif_info(split /;/, $motif_info_list)
        { 
            $count++;
            my $motif = "";
            my ($str_seq, $start, $end, $str_seq_length, $str_seq_motif_count) = ('', '', '', '', '');

            # 自定义区域信息
            if($motif_info =~ /^([ATCG]+)\((\d+)-(\d+)\)$/) 
            {
                $motif     = uc($1);
                my $aim_start = $2;
                my $aim_end   = $3;
                die "region must be in ascending order! $target_full_name\n" if($aim_start > $aim_end);
                ($str_seq, $start, $end, $str_seq_length, $str_seq_motif_count) = get_aim_str($target_full_name, $target_seq, $motif, $aim_start, $aim_end); # 获取指定目标区域motif的STR 
            }
            else
            {
                die "motif format info error, please read readme.txt : $target_full_name -> $motif_info\n";
            }
            
            # 检测是否重复
            my $title = "$motif,$start,$end";
            die "[Error] duplicate STR in $target_full_name: $motif $motif_info\n" if(exists($hashResult{$target_full_name}{$title})); # 重复STR
            
            # 信息保存
            $hashResult{$target_full_name}{$title}{'sort_order'}          = $count;      # 排序
            $hashResult{$target_full_name}{$title}{'original_info'}       = $motif_info; # 当前motif原信息
            $hashResult{$target_full_name}{$title}{'motif'}               = $motif; # motif
            $hashResult{$target_full_name}{$title}{'motif_length'}        = length($motif); # motif
            $hashResult{$target_full_name}{$title}{'str_seq_start'}       = $start; # str序列起始
            $hashResult{$target_full_name}{$title}{'str_seq_end'}         = $end;   # str序列终止
            $hashResult{$target_full_name}{$title}{'str_seq_length'}      = $str_seq_length;      # str序列长度
            $hashResult{$target_full_name}{$title}{'str_seq_motif_count'} = $str_seq_motif_count; # str序列含有的motif数量
            $hashResult{$target_full_name}{$title}{'target_seq_length'}   = $target_seq_length;   # 参考序列长度
            $hashResult{$target_full_name}{$title}{'str_mark'}            = "$motif($str_seq_motif_count)";   # str序列特殊表示形式
            $hashResult{$target_full_name}{$title}{'str_seq'}             = uc($str_seq);   #$motif x $str_seq_motif_count;   # str序列
            $hashResult{$target_full_name}{$title}{'target_seq'}          = $target_seq;   # 参考序列
            $hashResult{$target_full_name}{$title}{'target'}              = $target_full_name;   # 片段名称
        }
    }

    # (2) 输出到文件
    my @heads = ('target', 
        'motif', 
        'motif_length', 
        'str_mark',
        'str_seq_start', 
        'str_seq_end', 
        'str_seq_length', 
        'str_seq_motif_count', 
        'target_seq_length', 
        'str_seq',
        'target_seq',
        );
    
    open OUT, ">$output_dir/target.fasta.STR.txt";
    print OUT (join "\t", @heads) . "\n";
    foreach my $target_full_name(sort keys %hashResult)
    {
        foreach my $title(sort {$hashResult{$target_full_name}{$a}{'sort_order'} <=> $hashResult{$target_full_name}{$b}{'sort_order'}} keys %{$hashResult{$target_full_name}})
        {
            my @datas = map{ $hashResult{$target_full_name}{$title}{$_} } @heads;
            print OUT (join "\t", @datas) . "\n";
        }
    }
    close OUT;

    return "$output_dir/target.fasta.STR.txt";
}



# 创建目录
sub make_dir{
    my @dirs = @_;
    foreach my $dir(@dirs)
    {
        mkdir $dir if(not -e $dir);
    }
}

# 检验文件是否为空
sub is_file_ok{
    my @files = @_;
    my $isOK = 1;
    foreach my $file(@files)
    {
        $isOK = 0 if (not -e $file or not -f $file or -s $file == 0);
    }    
    return $isOK;
}

# 检验目录是否存在
sub is_dir_ok{
    my $dir = shift @_;
    my $isOK = 0;
    $isOK = 1 if(-e $dir and -d $dir);
    return $isOK;
}

# 提取fasta序列
sub read_fasta{
    my $fasta = shift @_;
    my %hashFasta;
    my $IN = Bio::SeqIO->new(-file => $fasta, -format=>'Fasta') or die "Could not open up file $fasta: $!";
    while(my $inSeq = $IN->next_seq)
    {
        my $id  = $inSeq->display_id;
        my $seq = $inSeq->seq;
        $hashFasta{$id}= $seq;
    }
    return %hashFasta;
}

# 获取指定目标区域motif的STR
sub get_aim_str{
    my $target_full_name = shift @_;
    my $target_seq = shift @_;
    my $motif      = shift @_;
    my $aim_start  = shift @_;
    my $aim_end    = shift @_;

    my $start          = $aim_start;
    my $end            = $aim_end;
    my $motif_length   = length($motif);
    my $str_seq        = substr($target_seq, $aim_start - 1, $aim_end - $aim_start + 1);
    my $str_seq_length = length($str_seq);

    my $str_seq_motif_count      = $str_seq_length / $motif_length; # STR序列中包含的motif数量
    # my ($target_seq_motif_count) = ($target_seq =~ s/$motif/$motif/ig); # 参考序列中包含的motif数量
    # my ($str_seq_motif_count_R)  = ($str_seq =~ s/$motif/$motif/ig); # STR序列中包含的motif数量, 根据实际获取的STR序列计算
    
    die "[Error] str length can not be divided by motif len!>$target_full_name\n$target_seq\n" if($str_seq_length % $motif_length!=0);# 非整数倍
    # die "[Error] aimed STR Region Error(motif count in region $aim_start-$aim_end not equal regionLength/motifLength):>$target_full_name\n$target_seq\n" if($str_seq_motif_count != $str_seq_motif_count_R);# 指定区域的motif数量跟实际不符
    # die "[Error] no motif exists in Seq :>$target_full_name\n$target_seq\n" if($target_seq_motif_count == 0 or $str_seq_motif_count_R == 0);# 参考序列中不含有motif
        
    return ($str_seq, $aim_start, $aim_end, $str_seq_length, $str_seq_motif_count);  
}


# blast+ 建库
sub build_blast_index{
    print "[process] build blast DB\n";

    my $target_fasta_db_dir = "$output_dir/target_fasta_db";
    my $target_fasta_db     = "$target_fasta_db_dir/seq.fa";
    make_dir($target_fasta_db_dir);
    system("cp $target_fasta $target_fasta_db");
    
    system("$SOFT_MAKEBLASTDB -in $target_fasta_db -dbtype nucl  -parse_seqids > $target_fasta_db_dir/log.txt 2>&1");
    return  $target_fasta_db;  
}


# 读取片段motif信息
sub read_motif{
    my $motif_file = shift @_;
    my $model      = (exists $_[0]) ? $_[0] : "full_name";
    
    die "Lost $motif_file\n" if(is_file_ok($motif_file) == 0);
    my %hashMotif;
    open MOTIF, $motif_file;
    my $line1 = <MOTIF>;
       $line1 =~ s/[\r\n]//g;
    my @heads = split /\t/, $line1;

    while(<MOTIF>)
    {
        $_=~s/[\r\n]//g;
        my @datas = split /\t/, $_;
        
        my %hashTmp = map{ ($heads[$_], $datas[$_])} (0..$#heads);

        my $target        = $hashTmp{'target'};
        my $motif         = $hashTmp{'motif'};
        my $str_seq_start = $hashTmp{'str_seq_start'};
        my $str_seq_end   = $hashTmp{'str_seq_end'};
        ($target) = split /[\|]/, $target if($model eq 'clean_name');

        my $aim_str       = "$motif,$str_seq_start,$str_seq_end";
        map{ $hashMotif{$target}{$aim_str}{$_} = $hashTmp{$_} } keys %hashTmp;              
    }
    close MOTIF;

    return %hashMotif;
}


# 获取序列中motif对应的最长STR
# sub get_seq_str{
    # my $target_full_name = shift @_;
    # my $target_seq = shift @_;
    # my $motif      = shift @_;
    # my $minrep     = shift @_; # 最小重复次数
    # my $type       = exists $_[0] ? $_[0] : 'Check';

    # my $motif_length = length($motif);
       # $minrep       = $minrep - 1;
    # my $regexp       = "(($motif)\\2{".$minrep.",})";

    # # 匹配最长str
    # my $bestlength = 0;
    # my $getinfo    = "";
    # while($target_seq =~ /$regexp/ig)
    # {
        # my $str_seq        = $1; 
        # my $str_seq_length = length($str_seq);
        # my $end            = pos($target_seq);
        # my $start          = $end - $str_seq_length + 1;
        # my $str_seq_motif_count = $str_seq_length / $motif_length; # STR序列中包含的motif数量
 
        # if($str_seq_length > $bestlength)
        # {
           # $bestlength = $str_seq_length;
           # $getinfo    = "$str_seq|$start|$end|$str_seq_length|$str_seq_motif_count";
        # }
        # pos($target_seq) = $start; # 重新设置检测起点，以单碱基步频移动  
    # }
    # my ($str_seq, $start, $end, $str_seq_length, $str_seq_motif_count) = ('', '', '', 0, 0);
       # ($str_seq, $start, $end, $str_seq_length, $str_seq_motif_count) = split /\|/, $getinfo if($getinfo =~ /\w/);

    # my ($target_seq_motif_count) = ($target_seq =~ s/$motif/$motif/ig); # 参考序列中包含的motif数量
    # die "no motif '$motif' exists in Seq >$target_full_name\n$target_seq \n" if($type eq 'Check' and ($target_seq_motif_count == 0 or $bestlength == 0) ); # 参考序列中不含有motif
    # my $str_seq_motif_perc = ($target_seq_motif_count == 0) ? 0 : sprintf "%0.4f", $str_seq_motif_count / $target_seq_motif_count;# str序列中motif占参考序列中motif数量的比例

    # return ($str_seq, $start, $end, $str_seq_length, $str_seq_motif_count, $str_seq_motif_perc);
# }
sub get_seq_str{
    my $target_full_name = shift @_;
    my $target_seq = shift @_;
    my $motif      = shift @_;
    my $minrep     = shift @_; # 最小重复次数
    my $type       = exists $_[0] ? $_[0] : 'Check';

    my $motif_length = length($motif);
       $minrep       = $minrep - 1;
    my $regexp       = "(($motif)\\2{".$minrep.",})";

    my $str_seq = $target_seq;
    # 获取起始位置（字符串的第一个字符的位置）
    my $start = 0;
    my $str_seq_length = length($str_seq);
    # 获取结束位置（字符串的最后一个字符的位置）
    my $end = $str_seq_length - 1;
    my $str_seq_motif_count = $str_seq_length / $motif_length; # STR序列中包含的motif数量
	
    my $str_seq_motif_perc = sprintf("%.4f", 100);

    return ($str_seq, $start, $end, $str_seq_length, $str_seq_motif_count, $str_seq_motif_perc);
}

sub get_sample_list{
    my $fastq_dir = shift @_;

    my @samples;
    opendir FASTQ_DIR, $fastq_dir;
    while(my $file = readdir FASTQ_DIR)
    {
        my ($sample) = $file =~ /(.*)_1.fq/;
        push @samples, $sample if(defined $sample);
    }
    close FASTQ_DIR;
    my $sample_list = join ",", @samples;
    return $sample_list;
}

# -----------------------------------
# 函数：在没有任何 overlap 的情况下，将R1和R2合并
# 参数： (1) R1.fastq，(2) R2.fastq，(3) 输出 merged.fastq
# -----------------------------------
sub merge_no_overlap_files {
    my ($fq1, $fq2, $merged_out) = @_;

    open my $fh1, "<", $fq1 or die "Cannot open $fq1: $!";
    open my $fh2, "<", $fq2 or die "Cannot open $fq2: $!";
    open my $out_fh, ">", $merged_out or die "Cannot write $merged_out: $!";

    while (1) {
        # 读取 R1 的 4 行和 R2 的 4 行
        my @r1_lines;
        my @r2_lines;
        for (1..4) {
            my $line1 = <$fh1>;
            my $line2 = <$fh2>;
            # 如果任一文件已读到末尾，就停止循环
            if (!defined $line1 or !defined $line2) {
                last;
            }
            chomp($line1);
            chomp($line2);
            push @r1_lines, $line1;
            push @r2_lines, $line2;
        }

        # 如果读出的行不足 4 条，说明文件结束
        last if (scalar(@r1_lines) < 4 || scalar(@r2_lines) < 4);

        # 分别取出各行信息
        my ($header1, $seq1, $plus1, $qual1) = @r1_lines;
        my ($header2, $seq2, $plus2, $qual2) = @r2_lines;

        # 对 R2 序列做反向互补（因为双端测序 R2 是反向读取）
        my $seq2_rc = reverse($seq2);
        $seq2_rc =~ tr/ACGTacgt/TGCAtgca/;
        # 对应质量值也要翻转
        my $qual2_rc = reverse($qual2);

        # 找重叠长度
        my $overlap_len = find_overlap($seq1, $seq2_rc);

        # 构造新 ID（去掉 @ 或 /1 /2 等）
        # 这里只是简单加了个 merged_ 前缀
        my $new_id = $header1;
        $new_id =~ s/^@//;     
        $new_id =~ s/\/[12]$//;     # 去掉结尾 /1 或 /2
        $new_id = "\@merged_$new_id"; 

        # # 根据 overlap 不同分情况处理
        # if ($overlap_len >= 10) {
            # # overlap>=10，FLASH已经合并，跳过不再输出。
            # next;
        # }
        # els
		if ($overlap_len == 0) {
            # 无重叠，直接拼接
            my $merged_seq  = $seq1 . $seq2_rc;
            my $merged_qual = $qual1 . $qual2_rc;

            print $out_fh "$new_id\n";
            print $out_fh "$merged_seq\n";
            print $out_fh "+\n";
            print $out_fh "$merged_qual\n";
        }
        else {
            #   0 < overlap_len < 10
            #   做“小重叠合并”。以下是最简单的做法：
            #   只取 seq1 + (seq2_rc去掉重叠部分) 作为合并序列
            #   重叠区在 seq1 中已包含，不重复写
            #   质量值也作同样处理

            # overlap_len 长度的序列，是 seq1 末端 & seq2_rc 的开头
            my $merged_seq  = $seq1 . substr($seq2_rc,  $overlap_len);
            my $merged_qual = $qual1 . substr($qual2_rc, $overlap_len);

            print $out_fh "$new_id\n";
            print $out_fh "$merged_seq\n";
            print $out_fh "+\n";
            print $out_fh "$merged_qual\n";
        }
    }

    close $fh1;
    close $fh2;
    close $out_fh;
}

# -----------------------------------
# 函数：检查 seq1 3' 和 seq2_rc 5' 的最大重叠
# 返回重叠长度，若无重叠则返回 0
# -----------------------------------
sub find_overlap {
    my ($s1, $s2) = @_;
    my $max_check = length($s1) < length($s2) ? length($s1) : length($s2);
    for (my $len = $max_check; $len > 0; $len--) {
        my $end_s1   = substr($s1, length($s1) - $len, $len);
        my $start_s2 = substr($s2, 0, $len);
        if ($end_s1 eq $start_s2) {
            return $len;
        }
    }
    return 0;
}