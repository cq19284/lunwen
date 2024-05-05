import pandas as pd
from scipy.stats import spearmanr

# 读取sRNA表达量数据
expr_df = pd.read_csv('rpm429.csv', index_col=0)

# 读取存活率数据
survival_df = pd.read_csv('new.csv', index_col=0)

# 确保expr_df的列名与survival_df的索引匹配，然后合并数据
# 注意：这里直接将survival_df的'存活率'列添加到expr_df的右侧，基于列名（样本ID）自动对齐
expr_df = expr_df.join(survival_df['SR'])
print(expr_df)
# 执行相关性分析
correlations = {}
for srna in expr_df.index[1:]:
    # 计算当前sRNA与所有样本的存活率之间的Spearman相关系数和P值
    corr, p_value = spearmanr(expr_df.loc[srna, 'SR'], expr_df.loc[srna, expr_df.columns[:-1]])
    correlations[srna] = {'Correlation': corr, 'P-value': p_value}

# 将字典转换为DataFrame以便查看和排序
correlation_results = pd.DataFrame.from_dict(correlations, orient='index')

# 筛选并打印出显著相关（例如，P值<0.05）的sRNA
significant_srnas = correlation_results[correlation_results['P-value'] < 0.05]
significant_srnas = significant_srnas.sort_values(by='Correlation', ascending=False)
print(significant_srnas)