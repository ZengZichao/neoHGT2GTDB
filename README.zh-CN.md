[English](README.md)（主）| **简体中文**


# neoHGT2GTDB: 🧬 基于 GTDB 参考框架的水平基因转移

> **这不是官方的 neoHGT / neoHGT2。** neoHGT2GTDB 是 [neoHGT2](https://github.com/cgneo/neoHGT) 的一个独立、非官方适配版本，在本仓库中维护。它未经上游作者或 [Qiyun lab](https://github.com/qiyunlab) 的认可，与其没有任何隶属关系，也不受其支持。官方项目仍位于 [cgneo/neoHGT](https://github.com/cgneo/neoHGT)；此处不使用任何官方 logo 或品牌标识。

**neoHGT2GTDB** 是对 [neoHGT2](https://github.com/cgneo/neoHGT) 的 GTDB 适配。它保留上游的 NCBI 路线、数据库下载器、搜索与分析流水线，并新增 GTDB 路线：由本地只读 [GTDB](https://gtdb.ecogenomic.org/) 发行版构建参考数据库、记录在只追加清单中的合成 taxId、GTDB 阶元与供体标签，以及可选的 NCBI+GTDB 混合树。GTDB 路线已在本地 R09-RS220 副本上完成端到端运行；官方平铺布局尚未通过完整建库验证（见 [GTDB 数据库（实验性）](#gtdb-数据库实验性)）。

命令行入口仍是 `neoHGT`，继承自上游，因此既有的 neoHGT/neoHGT2 脚本可以继续工作。GTDB 新增功能位于 `gtdb` 工具子命令和新的 `--source gtdb|hybrid` 数据库选项中。

**命名说明：** *仓库*和 *pip 发行包*名为 **neoHGT2GTDB**；*导入包*和 *命令*仍为 `neoHGT`。这是兼容性选择：文档示例和既有脚本都使用上游拼写。`pip install .` 会安装名为 `neoHGT2GTDB` 的发行包，它提供 `neoHGT` 命令和 `neoHGT` 包。

项目沿革：[HGTector](https://github.com/qiyunlab/HGTector) → neoHGT (v1) → [neoHGT2](https://github.com/cgneo/neoHGT) → **neoHGT2GTDB**（本仓库）。如果你有兴趣复现[论文](https://www.nature.com/articles/s41467-023-36487-z)中的结果，前一版本 neoHGT1 仍可在[这里](https://github.com/cgneo/neoHGT-v1)找到。继承自 neoHGT2 的 NCBI 路线问题最好报告到[上游 issue 跟踪器](https://github.com/cgneo/neoHGT/issues)；GTDB 适配中的问题应报告到本仓库的 issue 跟踪器，欢迎合作者加入。

## 目录

1. [关于](#关于)
2. [安装](#安装)
3. [GTDB 数据库（实验性）](#gtdb-数据库实验性)
4. [NCBI 数据库](#ncbi-数据库)
5. [编译](#编译)
6. [搜索](#搜索)
7. [分析](#分析)
8. [文档](#文档)
9. [卸载](#卸载)
10. [发表文献](#发表文献)
11. [文献引用](#文献引用)
12. [许可证](#许可证)

## 关于

**neoHGT2GTDB** 为感兴趣的物种预测水平/侧向基因转移（HGT）。与上游 [neoHGT](https://github.com/cgneo/neoHGT) 项目一样，它基于 Qiyun lab 开发的 [HGTector](https://github.com/qiyunlab/HGTector)。GTDB 路线在 GTDB 标准化、阶元完备的原核生物分类体系上重建参考框架，并可选择把 NCBI 树的顶端嫁接到 GTDB 树上。GTDB 路线与 NCBI 路线的统计行为正面对照尚未完成；GTDB 指南给出了现有实测数字和每种数据集的适用范围。

## 安装

本项目完全以 Python 3 编写。

第 1/3 步：

```sh
conda create -n neoHGT -c conda-forge python=3 pyyaml pandas matplotlib scikit-learn
conda activate neoHGT
```

第 2/3 步：从 bioconda 安装 DIAMOND。GTDB 路线已用 DIAMOND 2.1.10 测试。

```sh
conda install -c bioconda diamond
diamond version
```

第 3/3 步：安装本仓库。

```sh
conda activate neoHGT
pip install .
```

请在本仓库的克隆中运行。若要安装上游项目，请使用 `pip install git+https://github.com/cgneo/neoHGT.git`。

## GTDB 数据库（实验性）

GTDB 路线直接从**本地只读 GTDB 发行版**构建参考数据库。工具不会下载任何内容：分类体系由 GTDB 的 `*_taxonomy_r220.tsv` 和 `*_metadata_r220.tsv` 转换而来，蛋白序列从每个基因组的 `.faa` 文件读取，每个合成 taxId 都记录在只追加清单中，使 ID 在多次重建之间保持稳定。

### 环境

未配置工作区时，输出可以放在任何位置。设置 `NEOHGT_ROOT` 后，工具会拒绝写出该根目录之外；无论是否设置工作区，只读的 GTDB 输入都不会被写入。以下示例使用这些变量：

```sh
export GTDB_DIR=/path/to/GTDB_r220
export NEOHGT_ROOT=/path/to/workspace
export NEOHGT_CACHE="$NEOHGT_ROOT/cache"
export NEOHGT_TMP="$NEOHGT_ROOT/tmp"
```

也可以加载项目自带的封装；它根据脚本自身位置推导代码目录：

```sh
export GTDB_DIR=/path/to/GTDB_r220
. /path/to/neoHGT2GTDB/tools/env.sh
```

该封装还会读取 `tools/env.local`（被 git 忽略），用于机器特定路径。完整的环境与工作区说明见 [docs/GTDB-数据库.md](docs/GTDB-数据库.md)。

### 快速开始

```sh
# 0) 只读预检：识别布局并执行数据集自洽核对。
#    n_checks 随数据集变化；r220 参考副本为 19。
neoHGT gtdb probe --dir "$GTDB_DIR" \
  --out "$NEOHGT_CACHE/logs/layout_probe.json"

# 1) 建库：古菌每科一个代表基因组（约 564 个基因组）
neoHGT database --source gtdb --gtdb-dir "$GTDB_DIR" \
  --rank-sample family --domain archaea --seq-source local \
  --compile diamond --threads 16 \
  -o "$NEOHGT_CACHE/db/ar_family"

# 2) 搜索：DIAMOND 库已内嵌 GTDB taxId，staxids 直接可用
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

完整指南涵盖[数据准备](docs/GTDB-数据库.md#数据准备)、[probe](docs/GTDB-数据库.md#只读预检-probe)、[建库](docs/GTDB-数据库.md#建库-database)、[搜索](docs/GTDB-数据库.md#搜索-search)、[分析](docs/GTDB-数据库.md#分析-analyze)、[混合源](docs/GTDB-数据库.md#混合源-hybrid)、[资源预算](docs/GTDB-数据库.md#资源预算)和[故障排查](docs/GTDB-数据库.md#故障排查)。英文版见 [docs/en/GTDB-database.md](docs/en/GTDB-database.md)。

### 与 NCBI 路线不同的行为

- 远程分类信息获取被**强制关闭**（`--fetch-enable no`），因为合成 GTDB taxId 在 NCBI 没有记录。缺少 `staxids` 列的命中表会被拒绝（`--require-staxids` 默认 `yes`）；8 列标准 BLAST `outfmt 6` 表在 GTDB 模式下也会被拒收，因为它的 `qstart` 列会被读成 taxId。
- `nodes.dmp` 使用 NCBI 阶元词汇（`domain` 存为 `superkingdom`）并写满 13 列，这是 `diamond makedb` 的要求；命令行上 `domain` 和 `superkingdom` 都可接受。
- 供体以官方 GTDB 标签（`g__Halorubrum`）报告；`--donor-ncbi-map yes` 会写出 `hgts/<sample>.donor_ncbi.tsv`，离线地依据 metadata 的 `ncbi_taxid` 列把供体映射回 NCBI taxId。
- `neoHGT gtdb reconcile` 对 metadata、GTDB taxonomy 和本地 `.faa` 索引做差分，并写出五张 accession 级差异列表（即使为空也写）。任何差异都会以非零状态退出。
- `gtdb/` 下每个派生表都以 provenance 行开头：`# neoHGT version=…; gtdb_release=…; layout=…; generated=…`。`taxon.map.gz` 和 `genome.map.gz` 有意不带注释，因为 DIAMOND、BLAST 和 `read_prot2taxid` 会原样消费它们。
- GTDB-Tk 输出可以用 `neoHGT gtdb gtdbtk` 转换为库内合成 taxId，结果可直接用于 `neoHGT analyze --input-tax`。
- `neoHGT analyze --taxon-synonyms aliases.tsv`（`alias<TAB>current`）把 GTDB 旧名映射到当前名，作用于 `--self-tax`、`--close-tax`、`--tax-include`、`--tax-exclude` 和 `--input-tax`。
- `--keep-placeholder-names`（默认 `yes`）保留 `g__DUMJ01` 这类官方编号标签；r220 有 1,658 个类群使用这类标签。
- 在 GTDB 下，`-m precomp` 指向包含逐样本表的目录会被拒绝，因为上游循环会给每一批都喂第一张表。
- `--expand-members yes` 把代表基因组扩展为物种簇成员，并写出 `gtdb/members.tsv.gz`。
- `neoHGT gtdb tree-dist` 从发行版的逐阶元树计算 patristic 距离，并可用 `--db` 和 `--close-tax-out` 解析为合成 taxId。
- `--stream-extraction yes`（仅 NCBI 路线；默认关闭）用磁盘支持的实现替换上游蛋白抽取：实测每百万蛋白 882 → 138 MB，抽取时间为 2.3 倍，每百万蛋白约 100 MiB 临时盘。
- 以相同命令重新运行建库会复用 `db.faa` 和 `gtdb/extraction_state.json`，不会重复中断的抽取。
- `neoHGT search` 以**追加**方式写 `<output>/<sample>.tsv`；每次运行请使用新的 `-o` 目录，或先移走旧表。
- `db.dmnd` 随 taxId 的数值*范围*增长，因此 ID 基数默认 `100000000`（`--gtdb-taxid-base`，最小值 `10000000`）。
- `--source hybrid` 把 NCBI 树的顶端嫁接到 GTDB 树上。`--hybrid-ncbi-taxdump` 指向 NCBI 的 `nodes.dmp`/`names.dmp`（必需），`--hybrid-lca-policy root|domain|refuse` 决定两个来源都认领的蛋白如何处置。合并树在内存和 `gtdb/lineages.tsv` 中保留真实阶元，而 `nodes.dmp` 会把 DIAMOND 词表之外的阶元折叠为 `no rank`。可搜索的序列仍然只有 GTDB 蛋白。
- BLAST 路线**无法**用于 GTDB ID：`makeblastdb` 的本地 ID 上限为 50 字符，而 GTDB 序列 ID 为 151–270 字符。仅支持 DIAMOND。
- `--keep-faa drop` 会在建库成功后**删除 `db.faa`**。如果需要保留序列，请使用 `--keep-faa compress`。

### 不兼容的 CLI 缩写

argparse 的前缀匹配意味着新的长选项可能让旧缩写变得有歧义。以下缩写现在会报 “ambiguous option”，必须写出完整长选项：

| 子命令 | 现在失败的缩写 | 原先含义 |
| --- | --- | --- |
| `database` | `--a`、`--e`、`--geno`、`--l`、`--ra`/`--ran`、`--rep`、`--s`、`--ta`/`--tax`、`--ty`/`--typ`/`--type` | `--above`/`--latin`/`--genoids`/`--rank`/`--represent`/`--sample`/`--taxids`/`--typemater` 等 |
| `search` | `--r`/`--re`、`--s` | `--retries`、`--server` |
| `analyze` | `--t`/`--ta`/`--tax`、`--donor-n` | `--taxdump`、`--donor-name` |

`--h`、`--k`、`--w` 也是新的歧义缩写，但此前没有对应的长选项含义。`--g`/`--ge`/`--gen` 在上游就已在 `--genbank` 和 `--genoids` 之间有歧义，不是本次适配造成的。

### 实测开销

R09-RS220，全集含 113,104 个代表基因组 / 333,384,550 个蛋白：

| 层级 | 基因组数 | 蛋白数 | db.faa | db.dmnd | 墙钟时间 |
| --- | --- | --- | --- | --- | --- |
| 6 个基因组（冒烟） | 6 | 14,373 | 6.2 MB | 614 MB | ≈1 分钟 |
| 科级抽样，古菌 | 564 | 905,829 | 360 MB | 985 MB | 2 分 37 秒 |
| 一个科，`--reps` | 793 | 3,247,172 | 1.41 GiB | 2.07 GiB | 4 分 50 秒 |
| 属级抽样，细菌 | 23,112 | 约 6,000–8,000 万 | 28–37 GB | 42–110 GB | 外推 |
| 全部代表基因组集合 | 113,104 | 333,384,550 | 约 145 GiB | 约 1.5–3× db.faa | 抽取约 2–3 小时 + makedb 6–15 小时（外推） |

`db.faa` 每个蛋白占用 **466 字节**（三档实测 455 / 417 / 466），磁盘预检按此值估算，并为 `makedb` 留出 3× 峰值系数。

## NCBI 数据库

与 HGTector 相比，本项目继承的 neoHGT2 下载器修复了 HGTector 的一些已知问题，并使用多进程并行下载。

neoHGT2GTDB 保留上游的 `server.py`，它无需图形界面即可下载并安装数据库。

`server.py` 是**仅支持 NCBI** 的：它获取预构建的 `hgtdb_20230102` 数据库，对 GTDB 一无所知。GTDB（或 `--source hybrid`）数据库要通过命令行构建；参见上文的 [GTDB 数据库（实验性）](#gtdb-数据库实验性)。

推荐手动下载。

### Option1（手动）

#### （最全面、最新但耗时）

首先需要 `taxonmap`。有两种选择：

下面的链接会得到更小、更紧凑的 `taxdump`（解压后约 0.4 GB）。把链接粘贴到浏览器即可下载：

```text
https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.Z
```

下面的链接会得到最新、更全面的 `taxdump`（解压后约 1.8 GB）：

```text
https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.Z
```

___

接下来需要 `prot.accession2taxid`。有两种选择：

下面的链接会得到更小、更紧凑的 `prot.accession2taxid`（约 7 GB）：

```text
https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz
```

下面的链接会得到最新、更全面的 `prot.accession2taxid`（约 12 GB）：

```text
https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.FULL.gz
```

### Option2（手动）

#### （省时但不够全面、不够新）

[这个 OneDrive 链接](https://arizonastateu-my.sharepoint.com/personal/qzhu44_asurite_asu_edu/_layouts/15/onedrive.aspx?id=%2Fpersonal%2Fqzhu44%5Fasurite%5Fasu%5Fedu%2FDocuments%2FPublic%2FHGTector)上有一个预构建数据库。下载最新的 `hgtdb___.tar.xz` 文件后，可以跳过数据库下载与构建步骤，直接进入[编译](#编译)。

___

### 数据库下载

获得 `taxonmap` 和 `prot.accession2taxid` 后，运行 `neoHGT database` 下载约 291,000 个基因组文件并构建 DIAMOND 数据库（约 90 GB），具体取决于所选类别。

```sh
neoHGT database -c bacteria,archaea,viral -o <FOLDER_NAME>
```

`-c` 或 `--cats` 表示类别。下表说明自定义数据库的选项：

| 命令 | 选项 |
|---------|---------|
| -c (-cats) | archaea, bacteria, fungi, invertebrate, metagenomes, mitochondrion, plant, plasmid, plastid, protozoa, unknown, vertebrate_mammalian, vertebrate_other, 以及 viral（源自 [NCBI RefSeq 数据库的文件夹名](https://ftp.ncbi.nlm.nih.gov/genomes/refseq/)） |
| -o | 输出数据库目录 |

___

已知问题：如果看到下面的错误，不要慌张：

```text
ValueError: "___" is not a valid refseq genome category
```

出现这个错误，要么是类别名称有误（例如 NCBI 服务器上不存在 `virus` 文件夹，应使用 `viral`），要么是 NCBI 服务器暂时拒绝了连接。程序本身没有问题。

用以下命令测试连接：

```sh
rsync --list-only --no-motd rsync://ftp.ncbi.nlm.nih.gov/genomes/refseq/
```

如果看到 `rsync: failed to connect to ftp.ncbi.nlm.nih.gov: Connection refused (61)` 和 `rsync error: error in socket IO (code 10)`，说明存在网络问题，很可能在 NCBI 一侧。如果命令返回以 `drwxr-sr-x` 开头、随后是 `lrwxrwxrwx` 等条目的列表，说明连接已恢复，可以重新运行程序。

关于排除特定 taxId 等高级选项，请参考 HGTector 的[数据库](https://github.com/qiyunlab/HGTector/blob/master/doc/database.md)文档。

当程序以以下输出结束时：

```text
✅ Done
Combined protein sequences written to db.faa.
```

即可继续手动编译。

___

## 编译

编译只使用 DIAMOND，因为它经过优化、速度更快。

### 如果按照 Option1 操作

```sh
cd <DATABASE_DIRECTORY>
mkdir diamond
```

用你下载并解压好的文件替换标记为 `<‼️>` 的条目。确保这些文件都位于 `<DATABASE_DIRECTORY>` 的根目录下。下面为了清晰起见把命令拆成多行，实际使用时请保持为一行。

```sh
diamond makedb --threads 64 --in db.faa --taxonmap <‼️prot.accession2taxid.FULL> --taxonnodes <‼️taxdump>/nodes.dmp --taxonnames <‼️taxdump>/names.dmp --db diamond/db
```

### 如果按照 Option2 操作

第 1/2 步：

```sh
cd <DATABASE_DIRECTORY>
```

```sh
echo $'accession.version\ttaxid' | cat - <(zcat taxon.map.gz) > prot.accession2taxid.FULL
```

如果看到 `zcat: can't stat: taxon.map.gz (taxon.map.gz.Z): No such file or directory`，把 `taxon.map.gz` 重命名为 `taxon.map.gz.Z`。

第 2/2 步：

```sh
mkdir diamond
```

```sh
diamond makedb --threads 64 --in db.faa --taxonmap prot.accession2taxid.FULL --taxonnodes taxdump/nodes.dmp --taxonnames taxdump/names.dmp --db diamond/db
```

## 搜索

手动执行搜索被认为更安全，因为你将亲自运行 `diamond` 或 `blastp`。

第 1/3 步：

```sh
conda activate neoHGT
```

第 2/3 步：

```sh
diamond blastp --query <input_file.faa> --db <diamond_db> --outfmt 6 qseqid sseqid pident evalue bitscore qcovhsp staxids > <output_directory>/result.txt
```

第 3/3 步：

```sh
neoHGT search -i <input_file.faa> -m precomp -s <output_directory>/result.txt -t <taxdump_dir> -o <output_directory>
```

查看 `<output_directory>` 中的 `.tsv` 文件；如果非空，继续[分析](#分析)。

___

如果 `diamond` 无法工作，备选方案是 `blastp`，可在[这里](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html)找到。请按照[这篇文章](https://www.ncbi.nlm.nih.gov/books/NBK569841/)编译 BLAST 数据库。

```sh
makeblastdb -parse_seqids -in db.faa -out <blast_db_path>/<name_of_db> -title <name_of_db> -taxid_map <taxon.map>
```

`taxon.map` 文件可以通过以下命令获得：

```sh
gunzip -k taxon.map.gz
```

然后运行：

```sh
blastp -query <input_file.faa> -db <blast_db_path>/<name_of_db> -outfmt "6 qaccver saccver pident evalue bitscore qcovhsp staxids" > <output_directory>/result_blast.txt
```

最后：

```sh
neoHGT search -i <input_file.faa> -m precomp -s <output_directory>/result_blast.txt -t <taxdump_dir> -o <output_directory>
```

查看 `<output_directory>` 中的 `.tsv` 文件；如果非空，继续[分析](#分析)。

___

一站式搜索不推荐使用，因为 DIAMOND 依赖升级存在已知兼容性问题：

```sh
neoHGT search -i <input_file.faa.gz> -o . -m diamond -p 32 -d <diamond_db> -t <taxdump_dir>
```

详情见 [HGTector 的 issue 页面](https://github.com/qiyunlab/HGTector/issues/114)。

## 分析

遵循上游 neoHGT，图形输出格式为 `.pdf` 而非 `.png`，以获得更高分辨率并便于发表。

```sh
neoHGT analyze -i <output.tsv> -t <taxdump_dir> -o <output_dir>
```

这会生成若干文件，其中最重要的是 `scatter.pdf`，它是整条流水线的主要输出。

___

出于研究发表的严谨性，你可能希望完全掌控 `--self-tax` 和 `--close-tax`；否则程序会自动推断。

```sh
neoHGT analyze -i <output.tsv> -t <taxdump_dir> -o <output_dir> --self-tax <NCBI Taxonomy ID> --close-tax <NCBI Taxonomy ID>
```

___

## 文档

用户文档提供中英两版：

| 文档 | 语言 | 内容 |
| --- | --- | --- |
| [README.md](README.md) | 英文 | 本页：安装、GTDB 快速开始、NCBI 路线 |
| [README.zh-CN.md](README.zh-CN.md) | 中文 | 本页的中文镜像 |
| [docs/GTDB-数据库.md](docs/GTDB-数据库.md) | 中文 | 完整的 GTDB 指南 |
| [docs/en/GTDB-database.md](docs/en/GTDB-database.md) | 英文 | 完整的 GTDB 指南 |
| [docs/README.md](docs/README.md) | 英文 | 文档索引 |
| [docs/cli-reference.md](docs/cli-reference.md) | 英文 | 完整 CLI 选项参考 |
| [docs/config-reference.md](docs/config-reference.md) | 英文 | `config.yml` 参考 |
| [docs/troubleshooting.md](docs/troubleshooting.md) | 英文 | 常见问题与处理 |
| [CHANGELOG.md](CHANGELOG.md) | 英文 | 版本说明 |
| [CONTRIBUTING.md](CONTRIBUTING.md) | 英文 | 开发、测试与文档维护 |
| [CITATION.cff](CITATION.cff) | — | 机器可读的引用元数据 |

`docs/` 下的三个 JSON 文件（[`新增函数清单.json`](docs/新增函数清单.json)、[`方案编号清单.json`](docs/方案编号清单.json)、[`路由等价性清单.json`](docs/路由等价性清单.json)）是面向维护者的机器清单，保留中文文件名，没有英文版本；用途、字段和再生成命令见 [docs/README.md](docs/README.md)。

## 卸载

要卸载：

```sh
conda remove -n neoHGT --all
```

注意：`bioconda` 是 conda *通道*，不是包。要单独移除 DIAMOND，可运行 `conda remove -n neoHGT diamond`，或像上面那样删除整个环境。

## 发表文献

上游项目 neoHGT（v1）是作者在 LOB - Laboratoire d'Optique et Biosciences 实习期间的成果，该团队在 [Nature Communications](https://www.nature.com/articles/s41467-023-36487-z) 上发表了一篇论文；neoHGT2GTDB 将这项工作作为其科学基线和引用来源。

> Filée, J., Becker, H. F., Mellottee, L., Eddine, R. Z., Li, Z., Yin, W., Lambry, J.-C., Liebl, U., & Myllykallio, H. (2023, February 15). Bacterial origins of thymidylate metabolism in Asgard Archaea and Eukarya. Nature Communications. Retrieved February 18, 2023, from https://www.nature.com/articles/s41467-023-36487-z

## 文献引用

> Zhu Q, Kosoy M, Dittmar K. HGTector: [an automated method facilitating genome-wide discovery of putative horizontal gene transfers](https://bmcgenomics.biomedcentral.com/articles/10.1186/1471-2164-15-717). *BMC Genomics*. 2014. 15:717.

## 许可证

本仓库——非官方的 neoHGT2GTDB 适配——及其各上游项目均以 [BSD 3-clause](https://opensource.org/licenses/BSD-3-Clause) 许可证分发；参见完整的[声明](LICENSE)，其中载有每一条上游版权声明。

Copyright (c) 2026-, 曾子超 (Zeng Zichao, [ORCID 0000-0001-6553-970X](https://orcid.org/0000-0001-6553-970X)) 以及 neoHGT2GTDB 的各位贡献者。

Copyright (c) 2023-, Zhihui LI.

Copyright (c) 2013-, [Qiyun Zhu](mailto:qiyunzhu@gmail.com) 和 [Katharina Dittmar](mailto:katharinad@gmail.com)——上游 HGTector 和 neoHGT 的作者，本项目正是从这两个项目衍生而来。
