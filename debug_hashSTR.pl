#!/usr/bin/perl -w
use strict;
use warnings;
use Storable qw(store retrieve); # Import store and retrieve functions
use Data::Dumper; # Keep Data::Dumper for printing/debugging output

# 定义输入/输出目录
my $output_dir  = "/home/SSR/amp_bass/10ssr_216-634";  #  **修改为 output_dir 路径**
my @samples     = qw(FT100066721_L01_270-437-E01 FT100066721_L01_270-437-E02); # **修改为要调试的样本列表，保持和主脚本一致**

# File to store and load serialized %hashSTR (using Storable)
my $hashSTR_store_file = "$output_dir/hashSTR.dat"; # Changed filename to .dat extension
my $ploid;
$ploid = 2 unless defined $ploid;

# Load %hashSTR from file using Storable's retrieve
my %hashSTR;
if (-e $hashSTR_store_file) { # Check if the file exists before trying to load
    %hashSTR = %{ retrieve($hashSTR_store_file) }; # Retrieve and dereference
    print "Info: %hashSTR loaded from file: $hashSTR_store_file\n";
} else {
    die "Error: Serialized %hashSTR file not found: $hashSTR_store_file.  Make sure you run str_count_summary.pl first to create this file.\n";
}


print "Info: %hashSTR 从文件 '$hashSTR_store_file' 读取成功。\n";
# print Dumper(\%hashSTR); # 可选：打印加载后的 %hashSTR 内容，用于快速检查

# --------------------  从 str_count.pl 复制的后续代码 (动态调整 aim_str 和输出报告)  --------------------

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
                            # 如果满足条件 (超过 50%)，则执行非整数 motif_count 检测逻辑
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
    my $str_count_summary = "$output_dir/str_count-debug.txt";
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


# --------------------  调试代码结束  --------------------
print "Info: 调试脚本 debug_hashSTR.pl 执行完成。\n";