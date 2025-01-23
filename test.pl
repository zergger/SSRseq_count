#!/usr/bin/env perl
use strict;
use warnings;
use Bio::SearchIO;
use Data::Dumper; # 导入Data::Dumper模块

# 定义 BLAST 文件路径和目标查询名称
my $blast_file = '/home/SSR/amp_bass/10ssr_redo_656-778/mapping/FT100065359_L01_656-778-A01/FT100065359_L01_656-778-A01.blast';
my $target_query = 'FT100065359L1C012R00100537218:MERGE';

# 创建 Bio::SearchIO 对象
my $searchio = Bio::SearchIO->new(
    -format => 'blast',  # 根据实际 BLAST 格式选择，如 blast, blastxml 等
    -file   => $blast_file
);

# 遍历每一个查询
while ( my $result = $searchio->next_result ) {
    # 检查当前查询名称是否匹配目标查询名称
    if ( $result->query_name eq $target_query ) {
        print "查询名称: ", $result->query_name, "\n";
        print "查询长度: ", $result->query_length, "\n";

        # 遍历每一个命中（Hit）
        while ( my $hit = $result->next_hit ) {
            print "命中名称: ", $hit->name, "\n";
            print "命中长度: ", $hit->length, "\n";
            print "得分: ", $hit->score, "\n";

            # 遍历每一个比对（HSP）
            while ( my $hsp = $hit->next_hsp ) {
                print "  --- 比对区域 ---\n";
                print "  查询起始位置: ", $hsp->start('query'), "\n";
                print "  查询结束位置: ", $hsp->end('query'), "\n";
				my ($hit_start, $hit_end) = $hsp->range('hit') ;
                print "  命中起始位置: ", $hit_start, "\n";
                print "  命中结束位置: ", $hit_end, "\n";
                print "  比对长度: ", $hsp->length('total'), "\n";
                print "  比对得分: ", $hsp->score, "\n";
                print "  比对 E值: ", $hsp->evalue, "\n";
                print "  相似度: ", $hsp->percent_identity, "%\n";
				my $qerry_string = $hsp->query_string;
                print "  查询序列: ", $qerry_string, "\n";
				my $hit_string = $hsp->hit_string;
                print "  命中序列: ", $hit_string, "\n";
				print "  命中序列长度: ", length($hit_string), "\n\n";
                    my $get_seq     = ""; # 实际reads在目标STR区域的序列
                    my $before_seq  = ""; # 目标区域前面的序列
                    my $after_seq   = ""; # 目标区域后面的序列
                    my $target_pos  = $hit_start; # 参考序列上的绝对位置
my $str_seq_start = 65;
my $str_seq_end = 128;

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
print $gap_ref_pos, "\n";
        # if ($gap_ref_pos < $str_seq_start) {
            # # 缺口在目标区域之前，减少$str_seq_start
            # # 假设每个缺口减少1，可根据实际插入的碱基数量调整
            # $str_seq_start -= 1;
        # }
        # els
		if ($gap_ref_pos > $str_seq_start) {
            # 缺口在目标区域之后，增加$str_seq_end
            # 假设每个缺口增加1，可根据实际插入的碱基数量调整
            $str_seq_end += 1;
        }
print $str_seq_start, '====', $str_seq_end, "\n";
        # 如果缺口在目标区域内部，根据需求决定如何处理
        # 例如，可以选择将插入的碱基包括在$get_seq中，或者记录插入的位置
        # 此处不对内部缺口进行调整，直接处理
    }
    print "目标区域已根据缺口位置调整为：$str_seq_start 到 $str_seq_end\n";
} else {
    print "未检测到缺口，目标区域保持不变：$str_seq_start 到 $str_seq_end\n";
}
# 确保$str_seq_start和$str_seq_end的有效性
$str_seq_start = 1 if $str_seq_start < 1;
$str_seq_end = length($hit_string) if $str_seq_end > length($hit_string);
print Dumper(\@gap_positions);

                    for(my $hit_pos = 0; $hit_pos < length($hit_string); $hit_pos++)
                    {
                        my $hit_base   = substr($hit_string, $hit_pos, 1);
                        my $query_base = substr($qerry_string, $hit_pos, 1);
                        $get_seq    = $get_seq.$query_base    if($target_pos >= $str_seq_start and $target_pos <= $str_seq_end);
						print " get_seq: ", $get_seq, "\n";
                        $before_seq = $before_seq.$query_base if($target_pos < $str_seq_start);
                        $after_seq  = $after_seq.$query_base  if($target_pos > $str_seq_end);
                        # $target_pos++ if($hit_base=~/[ATCG]/i);# 只有当hit_base是碱基的时候，绝对位置向后移动一位，因为存在插入缺失，即：字符“-”
						$target_pos++;
                    }
            }
            print "\n";
        }
        # 如果只需要第一个匹配的查询，可以在此处退出循环
        last;
    }
}


