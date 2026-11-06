> 语言：简体中文 | English: [docs/en/GTDB-database.md](en/GTDB-database.md)
>
> README 中的精简版：[GTDB 数据库（实验性）](../README.zh-CN.md#gtdb-数据库实验性)

# GTDB 数据库参考

neoHGT2GTDB 可以从本地的只读 [GTDB](https://gtdb.ecogenomic.org/) 发行版构建参考数据库，并在该数据库上完成搜索与 HGT 分析。GTDB 提供标准化、层级完整的细菌与古菌分类，每个物种簇恰好有一个代表基因组；本工具把 GTDB 谱系转换为带稳定合成 taxId 的 NCBI 风格 taxdump，因此 DIAMOND、`neoHGT search` 与 `neoHGT analyze` 可以沿用原有流程。

GTDB 本身不含真核生物与病毒。需要把 NCBI 分类树的顶端（真核、病毒等）与 GTDB 原核树放在同一棵树中时，使用 `--source hybrid`；该模式需要另外准备 NCBI taxdump，并且当前版本只能搜索 GTDB 蛋白序列（见[混合源](#混合源-hybrid)）。

本页的性能数字均在本地 GTDB R09-RS220 副本上实测；标注为“外推”的条目是按实测吞吐和体积推算的，未实际执行。不同硬件、文件系统和数据集会得到不同结果。

## 目录

1. [概述](#概述)
2. [数据准备](#数据准备)
3. [环境与工作区](#环境与工作区)
4. [快速开始](#快速开始)
5. [只读预检 probe](#只读预检-probe)
6. [差异清单 reconcile](#差异清单-reconcile)
7. [构建 taxdump](#构建-taxdump)
8. [建库 database](#建库-database)
9. [搜索 search](#搜索-search)
10. [分析 analyze](#分析-analyze)
11. [其他工具](#其他工具)
12. [混合源 hybrid](#混合源-hybrid)
13. [输出文件与 meta.json](#输出文件与-metajson)
14. [断点续跑与 release 升级](#断点续跑与-release-升级)
15. [资源预算](#资源预算)
16. [已知限制](#已知限制)
17. [故障排查](#故障排查)
18. [附录：参数标定与实测细节](#附录参数标定与实测细节)

## 概述

GTDB 路线提供以下能力：

- 从本地只读 GTDB 发行版构建参考库，不需要联网下载。
- 把 GTDB 谱系转换为 NCBI 风格的 `nodes.dmp`/`names.dmp`，并为每个类群分配稳定的合成 taxId。
- 在搜索与分析阶段使用 GTDB 阶元（`domain`/`phylum`/`class`/`order`/`family`/`genus`/`species`）和官方标签（如 `g__Escherichia`）。
- 支持按代表基因组、类群子集、基因组清单或抽样方式选择建库基因组。
- 支持把代表基因组扩展为物种簇成员（`--expand-members`）。
- 支持把 GTDB-Tk 分类结果转换为库内合成 taxId，供 `analyze --input-tax` 使用。
- 支持使用 GTDB 的逐阶元树计算类群间 patristic 距离，并把结果作为 `analyze --close-tax` 的输入。
- 支持 `--source hybrid`，把 NCBI 分类树顶端嫁接到 GTDB 树上。

建库、搜索和分析是三个独立步骤。建库只需执行一次；搜索和分析可以对不同样本反复运行。

## 数据准备

### 发行版内容

建库需要以下本地文件：

| 内容 | 说明 |
| --- | --- |
| taxonomy 文件 | `*_taxonomy_rNNN.tsv` 或 `.tsv.gz`，每个 accession 一行 GTDB 谱系 |
| metadata 文件 | `*_metadata_rNNN.tsv`，包含代表基因组、质量分级、NCBI taxId 等列 |
| 蛋白序列 | 每个代表基因组一个 `.faa` 文件 |
| 可选：IDlist | `*_Sample_<rank>.IDlist`，用于快速选择类群子集 |
| 可选：树文件 | `tree_bac120/`、`tree_ar53/` 下的逐阶元折叠树 |

支持三种目录布局：

| 布局 | 特征 |
| --- | --- |
| `local-renamed` | 第三方整理过的副本：序列 ID 经过改写，文件分散在 `metadata_bac120/`、`metadata_ar53/` 等子目录 |
| `official-plain` | 官方 release 解包后的平铺布局：taxonomy 文件直接放在数据根目录 |
| `official-gz` | 与 `official-plain` 相同，但 taxonomy/metadata 为 `.gz` 压缩 |

`probe` 和 `taxdump` 会自动识别布局；也可以用 `--local-layout` 显式指定。`gtdb taxdump` 只需要 taxonomy 文件，`neoHGT database --source gtdb` 还要求 metadata 文件。

### 只读约定

`--gtdb-dir` 指向的数据集全程只读。工具在写入前检查输出路径，拒绝任何落在输入数据集内部的路径；构建结束后还会审计输入文件的 mtime。请不要把输出目录设在 GTDB 发行版内部。

### 显式指定文件

自动发现不适用于你的目录结构时，可以显式给出文件：

```sh
neoHGT gtdb probe --dir "$GTDB_DIR" \
  --gtdb-taxonomy "$GTDB_DIR/bac120_taxonomy_r220.tsv" \
  --gtdb-metadata "$GTDB_DIR/bac120_metadata_r220.tsv" \
  --out "$NEOHGT_CACHE/logs/layout_probe.json"
```

多个文件用逗号分隔。`--gtdb-release` 用于覆盖从文件名推断的 release 标签。

## 环境与工作区

### 两种使用方式

**不使用工作区（默认）**：输出可以写到任意有写权限的位置；工具只保护只读的 GTDB 输入数据集。适合已经用其他方式管理输出目录的用户。

**使用工作区（推荐）**：设置 `NEOHGT_ROOT` 后，所有输出必须位于该目录内。这样便于清理、迁移和审计，也能避免误写到数据集或系统目录。

### 环境变量

| 变量 | 用途 |
| --- | --- |
| `GTDB_DIR` | 本地 GTDB 发行版根目录（只读） |
| `NEOHGT_ROOT` | 工作区根目录；设置后所有输出必须位于其中 |
| `NEOHGT_CODE` | 源码目录；`tools/env.sh` 会自动推导 |
| `NEOHGT_CACHE` | 数据库、搜索结果、分析结果和日志的存放目录 |
| `NEOHGT_TMP` | 建库和搜索的临时目录，可安全清空 |
| `NEOHGT_PYTHON` | Python 解释器路径 |
| `NEOHGT_DIAMOND` | DIAMOND 可执行文件路径 |

`tools/env.sh` 是一个便捷封装：它根据脚本自身位置推导 `NEOHGT_CODE`，默认把工作区根设为代码目录的上一级，并创建缓存和临时目录。它不会写入工作区之外。

```sh
# 方式一：显式设置（适合脚本和集群任务）
export GTDB_DIR=/path/to/GTDB_r220
export NEOHGT_ROOT=/path/to/workspace
export NEOHGT_CACHE="$NEOHGT_ROOT/cache"
export NEOHGT_TMP="$NEOHGT_ROOT/tmp"
export DB="$NEOHGT_CACHE/db/ar_family"
export OUT="$NEOHGT_CACHE/analyze/sample"

# 方式二：使用项目自带的封装
export GTDB_DIR=/path/to/GTDB_r220
. /path/to/neoHGT2GTDB/tools/env.sh
```

如果不想在每个 shell 中重复设置，可以把 `GTDB_DIR` 和解释器路径写入仓库的 `tools/env.local`；该文件被 `.gitignore` 排除，不会提交。

本页后续示例使用 `$GTDB_DIR`、`$NEOHGT_CACHE` 和 `$NEOHGT_TMP`。把它们替换成你的实际路径即可；所有命令都可以直接复制执行。

## 快速开始

以下命令构成一条最小的 GTDB 建库、搜索和分析流程。

```sh
# 0) 只读预检：识别布局并执行数据集自洽核对
neoHGT gtdb probe --dir "$GTDB_DIR" \
  --out "$NEOHGT_CACHE/logs/layout_probe.json"

# 1) 建库：古菌每科一个代表基因组，约 564 个基因组
neoHGT database --source gtdb --gtdb-dir "$GTDB_DIR" \
  --rank-sample family --domain archaea --seq-source local \
  --compile diamond --threads 16 \
  -o "$NEOHGT_CACHE/db/ar_family"

# 2) 搜索：DIAMOND 库已内嵌 taxId，staxids 直接可用
diamond blastp --query sample.faa \
  --db "$NEOHGT_CACHE/db/ar_family/diamond/db" --threads 16 \
  --outfmt 6 qseqid sseqid pident evalue bitscore qcovhsp staxids \
  --max-target-seqs 25 --sensitive -o hits.tsv
neoHGT search -i sample.faa -m precomp -s hits.tsv \
  -t "$NEOHGT_CACHE/db/ar_family/taxdump" \
  -o "$NEOHGT_CACHE/search/sample"

# 3) 分析：使用 GTDB 阶元与官方供体标签
neoHGT analyze -i "$NEOHGT_CACHE/search/sample/sample.tsv" \
  -t "$NEOHGT_CACHE/db/ar_family/taxdump" \
  -o "$NEOHGT_CACHE/analyze/sample" \
  --self-rank genus --donor-rank genus --donor-label prefixed \
  --donor-ncbi-map yes
```

每一步的输入、输出和选项见后续章节。参考库的规模应覆盖多个科或更高级别；只含一个科的库缺少远缘类群，通常无法形成 distal 组。

## 只读预检 probe

```sh
neoHGT gtdb probe --dir "$GTDB_DIR" \
  --out "$NEOHGT_CACHE/logs/layout_probe.json"
```

`probe` 只读输入，输出 JSON 报告，内容包括：

- 布局判定、release、序列 ID 约定和 `_SerialNumber_` 锚点覆盖率；
- 基因组文件数、域集合、taxonomy/metadata 文件列表；
- 可用的 IDlist 档与代表基因组本地序列情况；
- 自洽核对结果 `selfcheck`，其中 `n_checks` 是本次数据集实际执行的检查数。

自洽核对的数量取决于数据集：每个 metadata 表、taxonomy 文件和 IDlist 档都会产生检查项。r220 参考副本为 19 项，但其他数据集可能不同；请以输出中的 `selfcheck.n_checks` 和 `selfcheck.status` 为准，不要把 19 当作固定值。

默认执行全部核对。跳过核对会降低可复现性保证：

| 子命令 | 跳过自检的选项 |
| --- | --- |
| `neoHGT gtdb probe` | `--no-selfcheck` |
| `neoHGT gtdb taxdump` | `--no-selfcheck` |
| `neoHGT database --source gtdb` | `--gtdb-no-selfcheck` |

`probe` 还有 `--no-deep`，用于跳过逐 accession 的对齐检查，适合在大型数据集上快速查看布局。输入不合格时命令以非零状态退出；`probe` 本身不写数据库。

## 差异清单 reconcile

`probe` 只报告计数是否相符。需要知道“哪些 accession 不符”时使用 `reconcile`：

```sh
neoHGT gtdb reconcile --dir "$GTDB_DIR" \
  --out "$NEOHGT_CACHE/reconcile" \
  --report "$NEOHGT_CACHE/reports/reconcile.json"
```

它在 metadata 全表、GTDB taxonomy 和本地 `.faa` 索引之间做集合差，输出五张 TSV：

| 文件 | 含义 |
| --- | --- |
| `missing_faa.tsv` | metadata 标记为代表、但本地没有序列文件的 accession |
| `extra_faa.tsv` | 本地有 `.faa`、但 metadata 不认识的 accession |
| `nonrep_faa.tsv` | 本地落盘的非代表（成员）基因组文件 |
| `metadata_not_in_taxonomy.tsv` | metadata 有、taxonomy 没有 |
| `taxonomy_not_in_metadata.tsv` | taxonomy 有、metadata 没有 |

每张表都有来源、release 和生成时间的注释头，即使为空也会写出。任意差异都会使命令以非零状态退出。r220 参考副本上五张表均为空；建库前建议对每个新数据集执行一次。

## 构建 taxdump

```sh
neoHGT gtdb taxdump --dir "$GTDB_DIR" \
  --out "$NEOHGT_CACHE/db/r220_taxdump"
```

产物为 `taxdump/nodes.dmp`、`taxdump/names.dmp` 和 `gtdb/` 下的映射表。r220 参考副本的实测结果：

- 146,380 个节点（GTDB 146,379 + root）；
- rank 计数：domain 2、phylum 217、class 607、order 2,030、family 5,460、genus 24,959、species 113,104；
- 转换约 42.9 秒，峰值内存约 1.5 GB，产物约 192 MB。

### 只有谱系文件的数据集

taxdump 是谱系串的函数，不需要 metadata 表。官方 release 的 taxonomy 文件可以直接放在数据根目录；`gtdb taxdump` 会识别这种布局并进入 `lineage_only` 模式：

```sh
neoHGT gtdb taxdump --dir /path/to/GTDB_r214/parent_directory \
  --gtdb-id-manifest /path/to/old_db/gtdb/taxid.map.gz \
  --out "$NEOHGT_CACHE/db/r214_taxdump"
```

该模式只执行不依赖 metadata 的核对（例如谱系深度恒为 7、同一 accession 不得有两个谱系），并打印被跳过的检查项。r214 参考副本实测约 113,318 个节点、7.9 秒。完整建库仍然要求 metadata 表。

谱系串不合格时，命令会在任何写入之前拒绝输入，并给出具体文件和行号。例如谱系深度不是 7、同一 accession 对应多个谱系，都会导致退出码 1，输出目录不会创建。

### 复用旧 ID 清单

升级 release 时，用 `--gtdb-id-manifest` 指向旧库的 `gtdb/taxid.map.gz`，可以保持已有类群的 taxId 不变。新 release 中消失的类群会写入 `gtdb/delnodes.txt`，其 taxId 保留在 `gtdb/reserved_ids.txt`，不会被后来的类群复用。`--gtdb-taxid-base` 必须与旧库一致（默认 `100000000`）。

## 建库 database

### 基本命令

```sh
neoHGT database --source gtdb --gtdb-dir "$GTDB_DIR" \
  --reps --seq-source local --compile diamond --threads 16 \
  -o "$NEOHGT_CACHE/db/r220_reps"
```

`--source` 在 `database` 阶段接受 `ncbi`、`gtdb`、`hybrid`，没有 `auto`。未给出时读取 `config.yml` 的 `database.source`（出厂值 `ncbi`），命令行始终优先。

### 常用选项

| 选项 | 作用 |
| --- | --- |
| `--gtdb-dir` | GTDB 发行版根目录；未给出时读 `$GTDB_DIR` 或 `config.yml` |
| `--gtdb-taxa` | 只包含这些 GTDB 类群的后代；接受 `f__Enterobacteriaceae` 或 GTDB taxId |
| `--reps` | 只使用物种簇代表基因组 |
| `--genome-list` | 从文件读取要包含的 GTDB accession，每行一个 |
| `--rank-sample <rank>` | 每个类群只取质量分最高的一个基因组 |
| `--sample N --rank <rank>` | 每个 `<rank>` 类群最多取 N 个基因组；可配合 `--above` |
| `--max-genomes N` | 截取选择结果的前 N 个基因组，便于冒烟测试 |
| `--expand-members yes` | 把代表基因组扩展为物种簇成员 |
| `--min-completeness` / `--max-contamination` | 按基因组完整度和污染度过滤 |
| `--mimag high\|medium\|low\|any` | 按 MIMAG 质量分级过滤 |
| `--type-material gtdb\|ncbi` | 按 type material 定义过滤；`--typemater` 是 `ncbi` 的别名 |
| `--domain bacteria\|archaea\|both` | 限制到一个域或同时构建两个域（默认 `both`） |
| `--seq-source local\|ncbi\|gtdb-tar\|ena` | 序列来源；本版本只有 `local` 有读取执行器 |
| `--compile none\|diamond\|blast\|both` | 编译哪种搜索库；GTDB 路线只支持 DIAMOND |
| `--threads N` | DIAMOND 编译线程数 |
| `--tmpdir DIR` | 临时目录；未给出时使用 `workspace.tmp_dir` 或系统临时目录 |
| `--keep-faa keep\|compress\|drop` | 建库后如何处理 `db.faa`，默认 `keep` |
| `--gtdb-taxid-base N` | 合成 taxId 起始值，默认 `100000000`，最小 `10000000` |
| `--gtdb-id-manifest PATH` | 复用旧库的 `gtdb/taxid.map.gz` |
| `--input-fingerprints yes` | 记录输入文件的 sha256 与大小 |
| `--workspace-out DIR` / `--workspace-tmp DIR` | 覆盖本次运行的输出/临时目录 |
| `--overwrite` | 覆盖已有产物并重新抽取 |
| `--force` | 跳过容量预检或强制执行不支持的组合 |

`--keep-faa drop` 会在建库成功后**删除** `db.faa`；如果还需要该文件，请使用 `keep` 或 `compress`。`compress` 会写出 `db.faa.gz` 并删除原文件。

默认值：`--seq-source` 默认 `local`；`--allow-network` 默认按 `no`（关闭）处理；`--gtdb-taxid-base` 默认 `100000000`。

### 抽样与 `--above`

`--sample N --rank <rank>` 先为每个 `<rank>` 类群取至多 N 个基因组。谱系没有到达该阶元的基因组不参与这一轮；给出 `--above` 后，程序从抽样阶元的上一级开始逐级补足每个粗类群到 N 个（到门级为止）。`--rank-sample <rank>` 是另一条低成本路径：每个类群只取质量分最高的一个，不使用 `--above`。

在 GTDB 代表集上，谱系通常完整到种，因此 `--above` 往往不会改变结果；它主要对谱系较浅的输入（GTDB-Tk 结果、成员基因组扩展）有意义。

### 质量类选项

r220 metadata 全表有 596,859 行（bacteria 584,382 / archaea 12,477），其中代表基因组 113,104 行。以下计数由生产选择代码在全表上统计：

| 选项 | 规则 | 全表 | 代表集 |
| --- | --- | ---: | ---: |
| `--type-material gtdb` | `gtdb_type_species_of_genus == 't'` | 4,545 | 3,347 |
| `--type-material ncbi` | `ncbi_type_material_designation` 非空且非 `na` | 22,256 | 16,371 |
| `--mimag high` | `mimag_high_quality == 't'` | — | 27,292 |
| `--mimag medium` | `mimag_medium_quality == 't'` | — | 85,752 |
| `--mimag low` | `mimag_low_quality == 't'` | — | 8 |
| `--mimag any` | 至少一个 MIMAG 分级为 `t` | 596,497 | 113,052 |

`--complete` 在 GTDB 路线会被拒绝：上游按 `ncbi_assembly_level` 判断“完整”，而该列描述的是序列提交方式，不是基因组完整度。请改用 `--min-completeness` 和 `--max-contamination`。r220 代表集上，completeness ≥ 50 保留 82,662 个，contamination ≤ 10 保留 82,464 个，两者同时满足保留 82,161 个。

`--mimag` 在缺少对应列的输入上会报错，而不是静默返回全表。

### 扩展到成员基因组

GTDB 每个物种簇只发布一个代表基因组。`--expand-members yes` 使用 metadata 的 `gtdb_genome_representative` 列把代表扩展为成员：

```sh
neoHGT database --source gtdb --gtdb-dir "$GTDB_DIR" \
  --gtdb-taxa f__Enterobacteriaceae --reps --expand-members yes \
  --compile none -o "$NEOHGT_CACHE/db/entero_members"
```

r220 参考副本上，793 个肠杆菌科代表扩展为 94,393 个成员基因组，写入 `gtdb/members.tsv.gz`。如果成员基因组没有本地蛋白文件，程序会在抽取前停止，并提示需要官方 `*_genomic_files_all` tar 或远程序列通道。

### 编译与容量

`--compile diamond` 使用 DIAMOND 编译 `diamond/db.dmnd`。GTDB 序列 ID 长度为 151–270 字符，超过 `makeblastdb` 的 50 字符上限，因此 `--compile blast` 或 `both` 默认被拒绝；`--force` 可以强制执行，但 BLAST 路线不可用。

建库前会做磁盘容量预检：空闲空间必须至少为估算 `db.faa` 的 3 倍，否则中止；`--force` 可以跳过预检。`--compile none` 只写出 `db.faa`、taxdump 和侧表，不生成 DIAMOND 库，`meta.json` 的 `diamond_version` 记为 `null`。

## 搜索 search

GTDB 库的 DIAMOND 数据库已经内嵌 taxId，因此命中表的 `staxids` 列可以直接使用。推荐“手工搜索 + 预计算结果”的路径：

```sh
diamond blastp --query q.faa --db "$DB/diamond/db" --threads 8 \
  --outfmt 6 qseqid sseqid pident evalue bitscore qcovhsp staxids \
  --max-target-seqs 25 --sensitive > hits.tsv

neoHGT search -i q.faa -m precomp -s hits.tsv -t "$DB/taxdump" \
  -o "$NEOHGT_CACHE/search/sample" --tmpdir "$NEOHGT_TMP/search/$$"
```

GTDB 模式的行为：

- 强制 `--fetch-enable no`，因为合成 taxId 在 NCBI 没有对应记录。
- `--require-staxids` 默认 `yes`；命中表缺少 `staxids` 列时直接报错。显式设为 `no` 时才降级为警告并继续。
- 8 列标准 BLAST `outfmt 6` 表会被拒绝：它没有 `staxids`，`qstart` 会被误读为 taxId。
- 无法确定 taxId 的命中写入 `unresolved_taxids.txt` 和 `unresolved_hits.tsv`。
- `search` 以追加方式写 `<output>/<sample>.tsv`。重复使用同一个 `-o` 目录会把每个蛋白写两遍；每次运行请使用新目录，或先移走旧表。

`--keep-placeholder-names` 默认 `yes`，来自 `config.yml` 的 `taxonomy.profiles.gtdb.keep_placeholder`。GTDB 有 1,658 个类群使用 accession 派生代码作为官方标签（如 `2-12-FULL-45-22`、`g__DUMJ01`、`s__MK2026 sp010841923`）。启用 `--tax-capital` 或 `--tax-latin` 时，这些标签会被豁免；设为 `no` 则取消豁免。

`--taxon-synonyms` 接受 `alias<TAB>current` 两列 TSV，用于把 GTDB 旧名（如 `p__Proteobacteria`）映射到当前名（`p__Pseudomonadota`）。别名表在解析 `--self-tax`、`--close-tax`、`--tax-include`、`--tax-exclude`、`--input-tax` 之前生效。

`--taxmap` 用于命中表没有 `staxids` 且序列 ID 无法解析的兜底路径。它需要明文、按第一列字节序排序的映射文件（通常是解压后的 `taxon.map.gz`），并用磁盘锚点索引，避免把全量映射读入内存。传入 `.gz` 会被拒绝。

## 分析 analyze

```sh
neoHGT analyze -i "$NEOHGT_CACHE/search/sample/sample.tsv" \
  -t "$DB/taxdump" -o "$NEOHGT_CACHE/analyze/sample" \
  --self-rank domain --donor-rank genus --donor-label prefixed \
  --donor-ncbi-map yes
```

### 阶元与标签

- `--self-rank`、`--donor-rank`、`--tax-unirank` 接受 `domain`，并自动映射到 taxdump 中存储的 `superkingdom`；也接受 `superkingdom`。
- `--donor-label prefixed` 是默认值，输出官方 GTDB 标签（如 `g__Escherichia`）。
- `--donor-ncbi-map yes` 写出 `hgts/<sample>.donor_ncbi.tsv`，用 metadata 的 `ncbi_taxid` 列把供体映射回 NCBI taxId，完全离线。

### 分组参数

`--source` 在 `search` 和 `analyze` 阶段接受 `auto`、`ncbi`、`gtdb`、`hybrid`，默认 `auto`。`auto` 从库目录检测来源；`gtdb` 找不到 GTDB 标记时直接报错；`hybrid` 找不到时打印警告并按上游 NCBI 路线运行。

GTDB 分组的默认 `closesize` 为 200（古菌档 `gtdb_archaea` 为 50），因为 GTDB 树比 NCBI 树密。域由每个蛋白最佳命中的域级祖先投票决定（严格多数，平局视为未判定）。日志中的 `GTDB profile: input domain = ...` 和 `Close group has N taxa.` 是实际生效值。

`--close-size` 默认按域取档：细菌 200，古菌 50。它是下限，不是配额。树上取不出足够多的类群时会被截断；跨样本比较必须看日志中的实际数量。`--bandwidth` 可以是 `auto`、`grid`、`silverman` 或 0.1–1.0 之间的数；`--noise` 的方向不保证单调，调参时应查看自己样本的曲线。

`--from-scores` 复用已写出的 `scores.tsv`，跳过分组阶段，并在日志中说明哪些分组参数本次不影响结果。两条路径对同一蛋白集合和排名应给出一致预测；浮点列可能因 `scores.tsv` 的 6 位有效数字在最后一位不同。

### close_size 实测

在 564 个古菌代表基因组、4,025 个查询蛋白的库上：

| close_size | 10 | 50 | 100 | 200 | 400 |
| --- | ---: | ---: | ---: | ---: | ---: |
| 预测 HGT 条数 | 371 | 371 | 132 | 0 | 0 |

在 200 和 400 档，close 组吞掉全部最佳命中，distal 组为空，日志给出 `WARNING: No hit is assigned to distal group.`。因此古菌档默认使用 50；更大的 close_size 需要接近全量代表集的库密度才有意义。

细菌样本的曲线并不一致。对 4,896 个“每科一个代表”的库：*Escherichia* 在 10/50/100/200 恒为 771，Alicyclobacillaceae 恒为 753，f__ATN3 恒为 714；Arcobacteraceae 在 10/50 为 436，到 100/200 直接归零，因为该属的 self 组只有 2 个类群，close_size 增大后近邻组跳到整个门并吞掉全部命中。因此“细菌档 200 可用”是逐样本性质，不是库的性质。

### 供体为 0

`donor=0` 表示远缘命中的 LCA 粗于 `--donor-rank` 请求的阶元，不是“没有供体”。在 E. coli 对 4,896 基因组跨科库的实测中，请求属级时 10.4% 可标注、科级 13.1%、域级 100%；`0` 的 LCA 节点中有 591 个在 superkingdom、65 个在 class、21 个在 phylum、14 个在 order。需要更高解析度时放宽 `--donor-rank`。

### 显式 close 组与输入分类

`--close-tax` 可以接受 `tree-dist --close-tax-out` 生成的 taxId 列表，实现系统发育定义的 close 组。`analyze` 会先把分类树裁剪到命中涉及的范围；完全不可达的 taxId 会被丢弃并告警。若全部不可达，则报错说明用法。

`--input-tax` 接受 `gtdb gtdbtk` 生成的 `genome_stem<TAB>taxId` 表，用于给查询基因组指定输入分类。

## 其他工具

### GTDB-Tk 结果转换

GTDB-Tk 不是运行时依赖；`gtdb gtdbtk` 只读取结果表：

```sh
neoHGT gtdb gtdbtk --in BINs.classify.max_summary \
  --db "$NEOHGT_CACHE/db/ar_family" \
  --out "$NEOHGT_CACHE/gtdbtk/input_tax.tsv" \
  --report "$NEOHGT_CACHE/gtdbtk/gtdbtk_report.json"
```

支持 v2 `classification`、v1 `classify_info`、de novo 的 `gtdbtk.lineage`/`lineage` 以及 7 列等级写法。taxId 取最深的可匹配前缀；库里没有的谱系判为 unmapped，不会凭空分配 ID。报告 JSON 包含 `matched_rank` 分布和每个未匹配基因组的原因。

### 供体清单映射 map-report

```sh
neoHGT gtdb map-report --db "$DB" --in hgt_list.txt \
  --out donor_ncbi.tsv
```

输入可以是 GTDB 标签或合成 taxId，每行一个。输出 TSV 给出 NCBI taxId 映射；映射不到的行记为 unmapped，不算失败。`--to` 默认 `ncbi`。

### 系统发育距离 tree-dist

```sh
neoHGT gtdb tree-dist --gtdb-dir "$GTDB_DIR" --rank genus \
  --focus g__QDMS2 --top-k 20 --out dist.tsv --report dist.json
```

也可以用 `--tree` 直接指定树文件。程序计算 patristic 距离，输出按距离排序的类群表和距离分布。`--gtdb-dir` 模式会依次尝试该阶元的两个域树，自动选中包含焦点类群的那一棵。

r220 属级树实测：

| 焦点 | 类群总数 | 中位距离 | ≤0.1 | ≤0.3 | ≤0.5 |
| --- | ---: | ---: | ---: | ---: | ---: |
| `g__QDMS2`（古菌） | 1,846 | 3.00 | 0 | 1 | 27 |
| `g__Halorubrum`（古菌） | 1,846 | 3.03 | 1 | 10 | 37 |
| `g__Escherichia`（细菌） | 23,111 | 2.36 | 33 | 111 | 158 |

把树距离作为 close 组：

```sh
neoHGT gtdb tree-dist --gtdb-dir "$GTDB_DIR" --rank genus \
  --focus g__QDMS2 --max-distance 0.5 --db "$DB" \
  --out close.tsv --close-tax-out close.txt
neoHGT analyze -i search/sample.tsv -t "$DB/taxdump" -o "$OUT" \
  --close-tax "$(cat close.txt)" --self-rank genus --donor-rank genus
```

`--db` 用库的 `gtdb/labels.tsv` 把标签解析为合成 taxId，并只保留库中真实存在的类群。在 564 基因组古菌库上，27 个 ≤0.5 的属中有 21 个在库内没有命中，会被告警丢弃；使用剩余类群作为 close 组后，预测从 371 条（9.2% 蛋白）提高到 849 条（21.1%），53 个 ar53 保守标记基因仍全部得 0 分。

## 混合源 hybrid

`--source hybrid` 把 NCBI 分类树的顶端嫁接到 GTDB 树上：

- GTDB 的 `d__Bacteria` 与 `d__Archaea` 改挂到 NCBI 的 cellular root（taxId 131567）。
- 按顶层 rank 保留若干真核与病毒顶层节点；其余 NCBI 节点不进入树。
- `--hybrid-ncbi-taxdump` 指向包含 `nodes.dmp`/`names.dmp` 的目录，必需。
- `--hybrid-lca-policy root|domain|refuse` 决定跨源共享蛋白的处置：
  - `root`：保留每个候选，LCA 落在 cellular root；
  - `domain`：优先落到该候选自己的域节点，多候选或跨域蛋白记入 `x_source_proteins`；
  - `refuse`：丢弃 taxId 不属于 GTDB 区间（`>= --gtdb-taxid-base`）的候选；没有可留候选时写 `na`。

```sh
neoHGT database --source hybrid --gtdb-dir "$GTDB_DIR" \
  --hybrid-ncbi-taxdump /path/to/ncbi_taxdump \
  --hybrid-lca-policy domain --reps --seq-source local \
  --compile diamond -o "$NEOHGT_CACHE/db/hybrid"
```

合并后的树保持单一根。r220 参考副本实测：纯 GTDB 树 146,380 个节点，嫁接 25 个 NCBI 节点后为 146,404 个节点；净增 24，因为其中一个嫁接节点就是根本身。两个域的父节点为 131567，细菌与真核的 LCA 落在 cellular root。

DIAMOND 会校验 `nodes.dmp` 的 rank 列。当前 NCBI taxdump 使用 `cellular root`、`acellular root`、`realm`、`clade` 等新式 rank，DIAMOND 不支持。`write_taxdump()` 会把这些 rank 在 `nodes.dmp` 中写成 `no rank`，而内存树、`gtdb/lineages.tsv` 和 `gtdb/labels.tsv` 保留真实 rank。被压平的节点记录在 `meta.json` 的 `taxonomy.hybrid.ranks_flattened_for_diamond` 和 `gtdb/build_report.md` 中。

hybrid 库当前只能搜索 GTDB 蛋白序列（`--seq-source local`）。嫁接进来的 NCBI 节点只提供树的语义，不贡献命中；要搜索真核/病毒蛋白，需要另外获取 NCBI 参考蛋白库。

## 输出文件与 meta.json

### 库目录结构

```text
<db>/db.faa
<db>/genomes.tsv
<db>/lineages.txt
<db>/taxon.map.gz
<db>/genome.map.gz
<db>/meta.json
<db>/taxdump/{nodes.dmp,names.dmp}
<db>/diamond/db.dmnd
<db>/gtdb/taxid.map.gz
<db>/gtdb/{labels.tsv,lineages.tsv,genome2taxid.tsv.gz,stem2taxid.tsv.gz,
           prefix2genome.tsv.gz,ncbi_map.tsv.gz,build_report.md,
           input_fingerprints.tsv,extraction_state.json,fetch_plan.tsv,
           members.tsv.gz,delnodes.txt,reserved_ids.txt}
```

`fetch_plan.tsv` 只在 `--seq-source ncbi|gtdb-tar|ena` 时出现；`members.tsv.gz` 只在 `--expand-members yes` 时出现；`delnodes.txt` 和 `reserved_ids.txt` 只在复用旧 ID 清单且确有节点变化时出现。

### meta.json

| 字段 | 内容 |
| --- | --- |
| `tool` / `source` / `seq_source` | 版本、数据源、序列通道 |
| `gtdb_release` / `local_layout` | release 标签与布局判定 |
| `id_convention` / `taxonomy` | 合成 taxId 规则、起始 base、各 rank 节点数、被取代节点数 |
| `extraction_resumed` | 本次是否复用了已有 `db.faa` |
| `input_provenance` | 输入根目录、只读审计信息、taxonomy/metadata 路径、序列目录、`input_sha256` |
| `counts` | 基因组/蛋白/去重/无 taxId/类群数/`db.faa` 体积/最长序列 ID 等 |
| `diamond_version` / `taxid_map_sha256` | DIAMOND 版本与 taxonmap 摘要 |
| `selection` | 本次选择方式与实际命中数 |
| `selfcheck` / `selfcheck_items` | 自洽核对状态与逐项结果 |
| `config_snapshot` / `build_started` / `build_finished` | 生效配置快照与时间戳 |

### provenance 头

`gtdb/` 下的派生表、`lineages.txt`、`input_fingerprints.tsv`、`fetch_plan.tsv`、`members.tsv.gz`、`delnodes.txt`、`reserved_ids.txt` 以及 reconcile 清单的第一行都是：

```text
# neoHGT version=2.0b3+gtdb1; gtdb_release=220; layout=local-renamed; generated=<timestamp>
```

版本号来自 `neoHGT.__version__`，会随安装版本变化。

`taxon.map.gz` 和 `genome.map.gz` 有意不带注释：DIAMOND、BLAST 和 `read_prot2taxid` 会原样消费它们。`gtdbtk` 生成的 `--input-tax` 文件同样不能带注释行。读取 `gtdb/` 内部表的函数会跳过 `#` 行。

## 断点续跑与 release 升级

### 断点续跑

抽取阶段会把计数器写入 `<db>/gtdb/extraction_state.json`。以相同命令再次运行到同一 `-o` 目录时：

- `db.faa`、状态文件和选择指纹都存在且一致：跳过抽取，并从 `taxon.map.gz`/`genome.map.gz` 重建衍生产物。
- 选择指纹不一致：不复用旧 `db.faa`，自动重新抽取。
- `db.faa` 存在但 `taxon.map.gz` 缺失：报错并提示用 `--overwrite` 重建。
- 需要强制重抽时使用 `--overwrite`。

### 升级到新 release

1. 对新数据集运行 `probe` 和 `reconcile`。`n_checks` 由数据集结构决定；如果 `idlist_vs_taxa[*]` 或 `rep_genomes_have_local_fasta` 失败，说明副本不完整，应停止。
2. 用旧库的 `gtdb/taxid.map.gz` 构建新 taxdump：
   ```sh
   neoHGT gtdb taxdump --dir /path/to/new_release \
     --gtdb-release 232 \
     --gtdb-id-manifest /path/to/old_db/gtdb/taxid.map.gz \
     --out "$NEOHGT_CACHE/db/r232_taxdump"
   ```
3. 正式建库时同样传入 `--gtdb-id-manifest`，并保持 `--gtdb-taxid-base` 一致。
4. 建库后检查 `meta.json` 的 `taxonomy.rank_counts` 和 `n_superseded`，抽查旧节点 taxId 未变。
5. 新 release 的标签可能变化；历史 `hgts/*.txt` 中的标签需要按新的 `gtdb/labels.tsv` 解释，`gtdb/ncbi_map.tsv.gz` 也应重新生成。
6. 库规模变化会影响 `close_size` 的可用区间；换库后应重新标定（见[附录](#附录参数标定与实测细节)）。

被删除类群的 taxId 会保留在 `reserved_ids.txt`，新类群不会复用。r214/r220 双向实测中，共有节点 ID 零改动、零冲突、零复用。

## 资源预算

| 档位 | 基因组 | 蛋白 | `db.faa` | `db.dmnd` | 墙钟时间 |
| --- | ---: | ---: | ---: | ---: | --- |
| 6 基因组冒烟 | 6 | 14,373 | 6.2 MB | 614 MB | 约 1 分钟 |
| 科级抽样，古菌 | 564 | 905,829 | 360 MB | 985 MB | 2 分 37 秒 |
| 一个科，`--reps` | 793 | 3,247,172 | 1.41 GiB | 2.07 GiB | 4 分 50 秒 |
| 属级抽样，细菌（外推） | 23,112 | 6,000–8,000 万 | 28–37 GB | 42–110 GB | 外推 |
| 全量代表集（外推） | 113,104 | 333,384,550 | 约 145 GiB | 约 1.5–3× `db.faa` | 抽取约 2–3 小时 + makedb 6–15 小时 |

`db.faa` 的体积基准为 **466 字节/蛋白**（三档实测 455 / 417 / 466）。容量预检按此值估算，并为 `makedb` 留出 3× 峰值系数。

DIAMOND 的两个约束：

1. `db.dmnd` 的体积随 taxId 的数值范围增长，而不是随节点数增长。同一份 13 MB 库在 `--gtdb-taxid-base 1e9` 下为 5.7 GB，在 `1e8` 下为 620 MB。因此默认使用 `100000000`，并拒绝低于 `10000000` 的取值。
2. `nodes.dmp` 的 rank 列必须使用 NCBI 词表，且必须写满 13 列；本实现会自动处理 GTDB 的 `domain → superkingdom` 映射和 NCBI 新式 rank。

跨卷读取带宽实测：单流 60 MB/s，8 并发 127 MB/s。建议分档构建（先古菌科级，再属级，最后全量），并在构建后用 `--keep-faa compress` 释放空间。

## 已知限制

- GTDB 只含细菌与古菌，不含真核和病毒。需要这些类群时使用 `--source hybrid`。
- hybrid 库当前只能搜索 GTDB 蛋白；NCBI 节点只提供树的语义。
- 远程序列通道（`--seq-source ncbi|gtdb-tar|ena`）只导出 `gtdb/fetch_plan.tsv` 供人工执行，本版本没有下载执行器。`--allow-network yes` 只是使用这些通道的前置确认，不会发起下载。
- 文件名内嵌谱系止于属级；种级以 metadata/taxonomy 为准。
- 本地 GTDB 副本可能经过第三方改名或 ID 重写；本工具会记录输入 sha256 并执行自洽核对，但没有官方 MD5SUM 可用于对照。
- `--keep-faa drop` 会删除 `db.faa`。如果需要保留序列，请使用 `compress`。
- GTDB taxId 是合成 ID，不能与 NCBI taxId 混用于同一棵历史结果树。跨 release 稳定性依赖 `gtdb/taxid.map.gz`。
- 官方平铺布局可以用于 `taxdump`，但完整建库仍要求 metadata 表；目前只在带 metadata 的本地副本上完成过完整建库。
- BLAST 路线不可用于 GTDB ID：`makeblastdb` 的本地 ID 上限为 50 字符，而 GTDB 序列 ID 为 151–270 字符。
- 与 NCBI 路线的正面对照尚未完成；本页的性能与分组数字来自 GTDB 库内部实测。
- 属级抽样和全量档的体积、时间为外推值，未实际执行。

## 故障排查

| 现象 | 原因 | 处理 |
| --- | --- | --- |
| `no metadata table ...` | 数据集只有 taxonomy 文件 | 建库需要 `--gtdb-metadata`；只建 taxdump 可以继续 |
| `does not have 7 ranks` | 谱系串被截断或用非标准分隔符 | 使用官方原始 taxonomy 文件 |
| `representative ... has no local FASTA` | 缺少代表基因组 `.faa` | 补齐序列文件，或改用 `--seq-source gtdb-tar\|ena` 生成抓取计划 |
| `Refusing to write ... outside workspace` | 设置了 `NEOHGT_ROOT`，输出路径不在其中 | 把输出移到工作区内，或调整 `NEOHGT_ROOT` |
| `Refusing to write ... input dataset is read-only` | 输出路径落在 `--gtdb-dir` 内 | 换一个输出目录 |
| `ambiguous option` | 长选项缩写与多个选项匹配 | 写出完整长选项 |
| `No hit is assigned to distal group` | close 组过大，吞掉全部命中 | 减小 `--close-size`，或扩大参考库 |
| `Invalid taxonomic rank` | NCBI taxdump 使用了新式 rank | 使用本工具生成的 taxdump；它会为 DIAMOND 压平不支持的 rank |
| `makeblastdb` 失败或 ID 被截断 | GTDB ID 超过 50 字符 | 使用 DIAMOND |
| `taxon.map.gz` 缺失 | 建库未完成或产物被移动 | 用 `--overwrite` 重建 |
| `db.faa` 不存在 | 使用了 `--keep-faa drop` | 重新建库，或从 `db.faa.gz` 解压（`compress` 档） |
| 内存不足 | 全量建库或 NCBI 抽取占用高 | 分档构建；NCBI 路线可用 `--stream-extraction yes` 降低峰值内存 |
| 磁盘不足 | 容量预检失败 | 释放空间，或显式 `--force`（不推荐） |

## 附录：参数标定与实测细节

### 流式抽取（NCBI 路线）

`--stream-extraction yes` 用磁盘支持的实现替换上游的蛋白抽取，仅适用于 NCBI 路线，默认关闭。实测 200k/1M 蛋白规模：

| 实现 | 200k 蛋白 | 1M 蛋白 | 边际成本 | 外推到 333.4M 蛋白 | 抽取 1M 蛋白耗时 | 临时盘 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 上游 | 365 MB | 1,071 MB | 约 882 MB/百万蛋白 | 约 287 GB | 20.3 秒 | 0 |
| 流式 | 193 MB | 304 MB | 约 138 MB/百万蛋白 | 约 45 GB | 45.9 秒 | 104 MiB |

内存降低约 6.4 倍，代价是约 2.3 倍抽取时间和每百万蛋白约 100 MiB 临时盘。`db.faa`、`genome.map.gz`、`taxon.map.gz` 与上游实现逐字节一致。

### `--taxmap` 磁盘索引

`--taxmap` 使用每 2048 行一个锚点的磁盘索引（约为数据量的 0.2%），查询时二分定位加局部扫描。实测 10,000 行、spacing 512 时生成 20 个锚点，索引文件 21 行；第二次打开复用缓存。映射文件必须是明文并按第一列字节序排序，`.gz` 会被拒绝。

### 参数敏感性

- `--bandwidth` 是影响最大的旋钮之一，幅度依赖样本。halo 样本的 `auto`/`grid`/`silverman` 为 371/145/222 条预测；bathy 为 367/346/366。
- `--noise` 的方向不保证单调。halo 在 25/50/75 下为 374/371/325，bathy 在 25/75 下为 285/421；在双峰或长尾分布上方向可能翻转。
- `--input-cov` 主要影响输入分类，而不是预测条数。命中分类干净时 50–95 通常平坦；覆盖率接近阈值时，提高阈值会把输入分类从种退到属。
- `--from-scores` 与完整分析在预测集合和排名上一致；浮点列可能因 `scores.tsv` 的精度在最后一位不同。

### 测试夹具

`neoHGT/tests/data/GTDB/` 是仓库内的 KB 级夹具（约 118 KiB），包含 6 个基因组的蛋白摘录、30 条 `db.faa` 记录、32 个 taxdump 节点和 168 行命中样例。它可以直接作为一个只读的可搜索库使用，相关离线端到端测试见 `CONTRIBUTING.md`。
