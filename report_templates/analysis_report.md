# scAgent 单细胞分析报告

- 日期: {date}
- 任务: {query}
- 物种: {species} ｜ 平台: {platform} ｜ 组织: {tissue} ｜ 样本数: {n_samples}

## 分析路线

{narrative}

- Skills: {skills}
- 步骤: {route}

## 质控策略

{qc_protocol}

## 注释策略

{annotation}

## 审查

- 通过: {review_passed}
- 问题: {review_issues}

## 图

{figures}

每张图需要图注：画了什么、参数、能得出什么、不能得出什么。未执行的分析不得写成结果。

## 运行

- 执行成功: {execution_ok}

```
{stderr}
```

## 局限

统计推断以生物学重复为单位。UMAP 是可视化。注释是分层证据而非单一标签。
