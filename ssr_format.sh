# 假设 $1 是输入文件名（str_count.txt），$2 是输出文件名（processed_str_count.txt）
awk '
BEGIN { FS=OFS="\t" }
NR == 1 {
    # 输出表头，并删除 AimSTR 列
    $2 = "";
    print $0
}
NR > 1 {
    # 提取 Target 列中的数字部分和特征
    split($1, parts, "\\|");
    id = parts[1];
    features = parts[2];

    # 分割特征并统计特征种类
    split(features, feature_list, ";");
    feature_count = 0;
    delete unique_features;
    for (i in feature_list) {
        feature = gensub(/[\(\)-]/, ",", "g", feature_list[i]);  # 去掉括号及连字符
		gsub(/^[ \t]+|[ \t]+$/, "", feature);  # 去除首尾空格
		gsub(/,+$/, "", feature);  # 去除末尾的一个或多个逗号
        if (!(feature in unique_features)) {
            unique_features[feature] = ++feature_count;
        }
    }

    # # 打印 unique_features 数组以进行调试
    # for (key in unique_features) {
        # print "unique_features[" key "] = " unique_features[key] > "/dev/stderr"
    # }

    # 如果特征种类超过1，则根据 AimSTR 列确定新ID
    if (feature_count > 1) {
        feature = $2;  # 根据 AimSTR 列确定特征
		gsub(/^[ \t]+|[ \t]+$/, "", feature);  # 去除首尾空格
		# print "Current feature: " feature > "/dev/stderr"  # 打印当前特征
        new_id = id "_" unique_features[feature];
    } else {
        new_id = id;
    }

    # 替换 Target 列
    $1 = new_id;

    # 删除 AimSTR 列
    $2 = "";

    # 输出修改后的行
    print
}' "$1" > intermediate.txt

# 删除每行前3列中的空列
awk '
BEGIN { FS=OFS="\t" }
{
    for (i=1; i<=3; i++) {
        if ($i == "") {
            for (j=i; j<NF; j++) {
                $j = $(j+1)
            }
            NF--
            i--
        }
    }
    print
}' intermediate.txt > "$2"

# 删除中间文件
rm intermediate.txt

echo "Processing complete. Output saved to $2."
