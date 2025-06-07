# Miyabi システム利用ガイド

## システム概要

Miyabi は東京大学情報基盤センターと筑波大学計算科学研究センターが共同運営するスーパーコンピューターシステムである。以下の 2 つの主要なシステムで構成される：

- **Miyabi-G**: NVIDIA GH200 Grace-Hopper Superchip 搭載の演算加速ノード
- **Miyabi-C**: Intel Xeon Max 9480 を搭載した汎用 CPU ノード（FUJITSU Server PRIMERGY CX2550 M7）

## 利用の流れ

### 1. システムへのアクセス

- 接続先: `miyabi-c.jcahpc.jp`
- アクセス方法: SSH/SCP（2 要素認証必須）
- 認証方式: 公開鍵認証 + ワンタイムパスワード（OTP）

### 2. ファイル転送方法

以下の方法でファイル転送が可能：

- UNIX/Mac/WSL 環境: `scp`または`sftp`コマンド
- Windows 環境: WinSCP または WSL 経由の`rsync`

### 3. 実行環境

- ログインノード: プログラムの作成・編集、コンパイル、ジョブ投入用

  - メモリ制限: ユーザーあたり 16GB どの程度かはわからない。
  - 高負荷処理は禁止（プリポストノードを使用）どの程度か...以下略

- ファイルシステム:
  - ホーム領域: `/home`
  - 作業領域: `/work`

### 4. プログラミング環境

Intel oneAPI (2023.2.0) が利用可能：

| 種別 | 言語    | コマンド | 最適化オプション              | OpenMP   |
| ---- | ------- | -------- | ----------------------------- | -------- |
| 逐次 | Fortran | ifort    | -axSAPPHIRERAPIDS,CORE-AVX512 | -qopenmp |
|      | C       | icc      | -axSAPPHIRERAPIDS,CORE-AVX512 | -qopenmp |
|      | C++     | icpc     | -axSAPPHIRERAPIDS,CORE-AVX512 | -qopenmp |
| MPI  | Fortran | mpiifort | -axSAPPHIRERAPIDS,CORE-AVX512 | -qopenmp |
|      | C       | mpiicc   | -axSAPPHIRERAPIDS,CORE-AVX512 | -qopenmp |
|      | C++     | mpiicpc  | -axSAPPHIRERAPIDS,CORE-AVX512 | -qopenmp |

### 5. ジョブ実行

#### ジョブ種別

- バッチジョブ: スクリプトによる一括実行 こっちを使うのかな？
- インタラクティブジョブ: 対話的実行

#### 主要なキュー

- バッチジョブ用: `debug-c`, `short-c`, `small-c`など
- インタラクティブジョブ用: `intract-c_n1`, `intract-c_n2`
- プリポスト処理用: `prepost`, `prepost1_n1`

#### ジョブ管理コマンド

- `qsub`: ジョブ投入
- `qstat`: ジョブ状態確認
- `qdel`: ジョブ削除
- `qhold`: ジョブ保留
- `qrls`: 保留解除

### 6. バッチジョブスクリプト例

#### 単一ノード逐次実行

```bash
#!/bin/sh
#PBS -q debug-c
#PBS -l select=1
#PBS -l walltime=30:00
#PBS -W group_list=group1
#PBS -j oe
cd ${PBS_O_WORKDIR}
./a.out
```

#### MPI 並列実行（2 ノード、ノードあたり 48 プロセス）

```bash
#!/bin/sh
#PBS -q regular-c
#PBS -l select=2:mpiprocs=48
#PBS -l walltime=12:00:00
#PBS -W group_list=group1
#PBS -j oe
cd ${PBS_O_WORKDIR}
mpiexec.hydra ./a.out
```

## その他の機能

- **Open OnDemand**: ブラウザベースの利用環境

  - 仮想デスクトップ
  - JupyterLab
  - ジョブ管理機能

- **コンテナ環境**: Singularity による実行環境の提供
