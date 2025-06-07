Miyabiシステムのご利用に関心をお寄せいただきありがとうございます。初回登録（アカウント、SSH鍵登録、OTP設定など）が完了済みであることを前提として、Miyabi-Cシステムにおけるファイルのアップロードからジョブ実行までの手順について、ソースに基づき詳細にご説明します。

Miyabiシステムは、東京大学情報基盤センターと筑波大学計算科学研究センターが共同運営するスーパーコンピューターシステムです。本利用手引書は、その利用方法を説明した資料です。

Miyabiシステムは、NVIDIA GH200 Grace-Hopper Superchipを搭載したMiyabi-G (演算加速ノード) と、Intel Xeon Max 9480を搭載したMiyabi-C (汎用CPUノード) で構成される超並列クラスタ型スーパーコンピュータです。Miyabi-CはFUJITSU Server PRIMERGY CX2550 M7で構成され、Intel Xeon MAX 9480 CPUを2基搭載しています。

### 1. ファイルのアップロード (ファイル転送)

Miyabiシステムへのファイル転送には、SSHアクセスとしてssh/scpが利用可能です。外部から直接アクセスできるノードは、Miyabi-Cログインノードなどです。Miyabi-Cログインノードの接続先アドレスは **miyabi-c.jcahpc.jp** です。

ファイル転送には、以下の方法があります。

*   **scpまたはsftpコマンド**：UNIX/Mac環境やWSL (Windows Subsystem for Linux) 環境から標準的なコマンドでファイル転送が可能です。
    *   例 (scpコマンド): `scp /path/to/local/file username@miyabi-c.jcahpc.jp:/path/to/remote/directory`
*   **WinSCPを利用する場合 (Windows環境)**：WinSCPは、ファイルをドラッグ＆ドロップすることでファイル転送できるGUIクライアントです。WinSCPを使用する前に、SSH初回ログイン時にOTP認証を設定する必要があります。ログイン時にホスト名に `miyabi-c.jcahpc.jp`、ポート番号に `22`、ユーザー名、パスフレーズを入力し、ワンタイムパスワードを入力します。
*   **rsyncコマンドを利用する場合 (Windows環境)**：WSLを使用してLinux環境を構築することで、rsyncコマンドによるファイル転送も可能です。

### 2. ログイン環境とソフトウェアの利用

Miyabi-Cシステムを利用するには、まずMiyabi-CログインノードにSSHでログインします。Miyabi-Cログインノードはx86_64アーキテクチャです。SSHログインには、公開鍵認証とワンタイムパスワード（OTP）認証を併用する2要素認証を採用しています。対応しているターミナルソフトが必要です。初回SSHログイン時にOTP認証の設定が必要です。

ログインノード上では、プログラムの作成・編集、コンパイル/リンク、バッチジョブの操作、インタラクティブジョブの実行などが行えます。ただし、ログインノードの資源は多くのユーザーで共有するため、**高負荷な処理は行わないでください**。高負荷な前処理・後処理はプリポストノードを利用してください。Miyabi-Cログインノードで使用できるメモリ量は、ユーザーあたり16GBです。

Miyabiシステムでは、全ての計算ノードとログインノードにて、**共有ファイルシステム**を共有します。共有ファイルシステムには、ホーム領域(/home)と作業領域(/work)があります。Miyabi-C計算ノードからは、Ipomoea-01システムのストレージへはアクセスできません。

Miyabi-Cで利用可能なソフトウェアは、moduleコマンドやshow_moduleコマンドで確認できます。`module avail`で設定可能なアプリケーションやライブラリの一覧を表示できます。`module load <module>`で必要な環境変数を設定します。

### 3. プログラムのコンパイル/リンク (Miyabi-C向け)

Miyabi-Cログインノードでは、Miyabi-C計算ノード向けに最適化されたロードモジュールを作成することが可能です。Miyabi-CではIntel oneAPI (intel/2023.2.0) が利用可能です。

Intel oneAPIのコンパイル/リンクコマンドは以下の通りです。

| 種別       | 言語処理系 | コマンド     | 最適化オプション注1              | OpenMP注2 |
| :--------- | :--------- | :----------- | :------------------------------- | :-------- |
| 逐次(非MPI) | Fortran    | ifort        | -axSAPPHIRERAPIDS,CORE-AVX512    | -qopenmp  |
|            | C          | icc          | -axSAPPHIRERAPIDS,CORE-AVX512    | -qopenmp  |
|            | C++        | icpc         | -axSAPPHIRERAPIDS,CORE-AVX512    | -qopenmp  |
| MPI並列注3  | Fortran    | mpiifort     | -axSAPPHIRERAPIDS,CORE-AVX512    | -qopenmp  |
|            | C          | mpiicc       | -axSAPPHIRERAPIDS,CORE-AVX512    | -qopenmp  |
|            | C++        | mpiicpc      | -axSAPPHIRERAPIDS,CORE-AVX512    | -qopenmp  |

*   注1: Miyabi-C計算ノード向けに最適化されたロードモジュールを作成するための推奨オプションです。
*   注2: OpenMPオプションはデフォルトでは無効です。
*   注3: MPI並列を使用する場合はIntel MPIライブラリのmoduleのロードが必要です。

**コンパイル/リンクの例**:

*   **逐次プログラム (C)**: `icc -axSAPPHIRERAPIDS,CORE-AVX512 seq.c`
*   **スレッド並列プログラム (C, OpenMP)**: `icc -axSAPPHIRERAPIDS,CORE-AVX512 -qopenmp omp.c`
*   **MPI並列プログラム (C)**: `mpiicc -axSAPPHIRERAPIDS,CORE-AVX512 mpi.c`
*   **ハイブリッド並列 (C, MPI+OpenMP)**: `mpiicc -axSAPPHIRERAPIDS,CORE-AVX512 -qopenmp mpi_omp.c`

Intel oneAPIでは、LLVMベースの新しいコンパイラ(ifx/icx/icpx)も利用可能です。

FFTW, GSL, SuperLU, METIS/ParMETIS, Scotch/PT-Scotch, HDF5, Parallel HDF5, netCDF, Parallel netCDF, PETScなどのライブラリも利用可能です。これらのライブラリを利用するには、通常、事前にmoduleコマンドで環境を設定し、コンパイル/リンク時に適切なオプションを指定する必要があります。Miyabi-C計算ノード向けの例が各ライブラリのセクションに記載されています。

### 4. ジョブ管理システムと実行

プログラムの実行は、ジョブ管理システムによりバッチジョブまたはインタラクティブジョブとして処理を依頼します。

**ジョブの種類**:

*   **バッチジョブ**: スクリプト単位でジョブを実行。ノードダウン時などに自動再実行が可能です（オプションによる）。
*   **インタラクティブジョブ**: ユーザー端末からの入力により会話形式でジョブを実行。主にデバッグ等での利用を想定。ノードダウン時に再実行されません。

Miyabiシステムで利用可能なジョブ操作コマンドは以下の通りです。

*   `qsub`: ジョブ投入
*   `qstat`: ジョブ参照
*   `qdel`: ジョブ削除
*   `qhold`: ジョブ保留
*   `qrls`: ジョブ保留の解除

**Miyabi-Cで利用可能なキュー**:

*   **バッチジョブ用キュー**: ノード専有利用のキュー (`debug-c`, `short-c`, `small-c` など) があります。ルーティングキュー `regular-c` を指定すると、ノード数に応じて自動的にキューが選択されます。
*   **インタラクティブジョブ用キュー**: (`intract-c_n1`, `intract-c_n2`)。ジョブを実行してもトークンは消費されません。1ユーザーあたり1本のジョブのみ投入・実行可能です。
*   **プリポスト処理用キュー**: (`prepost`, `prepost1_n1` など)。ログインノードでの高負荷な作業に使用します。インタラクティブジョブとして利用可能で、トークンは消費されません。ログインノード(Miyabi-C)と同一環境です。1ユーザーあたり1本のみジョブの投入・実行可能です。

**主なジョブ投入オプション (-qsub コマンド)**:

*   `-q <name>`: キュー名/ルーティングキュー名を指定します。
*   `-l select=<num_nodes>:mpiprocs=<num_procs_per_node>:ompthreads=<num_threads>`: ジョブが利用するノード数、ノードあたりのMPIプロセス数、1プロセスあたりのスレッド数などを指定します。
*   `-l walltime=<HH:MM:SS>`: ジョブの実行時間制限を指定します。
*   `-W group_list=<projectcode>`: ジョブ実行時にトークンを消費するプロジェクトを指定します（必須）。
*   `-o <filename>`: 標準出力を指定ファイルに出力します。同名ファイルが存在する場合は追記モードで出力されます。
*   `-e <filename>`: 標準エラー出力を指定ファイルに出力します。同名ファイルが存在する場合は追記モードで出力されます。
*   `-j oe`: 標準エラー出力を標準出力にマージして出力します。
*   `-N <jobname>`: ジョブ名を指定します。
*   `-I`: インタラクティブジョブとしてジョブを投入します。

Miyabi-Cへのバッチジョブ実行時には、自動でIntel OneAPI (intel および impi) がloadされるため、バッチジョブスクリプト内で `module load intel`, `module load impi` を記述する必要はありません。

**Miyabi-C向けバッチジョブスクリプト例**:

プログラムファイル `a.out` を実行する場合の例です。

*   **単一ノード逐次ジョブ**:
    ```bash
    #!/bin/sh
    #------ qsub option --------#
    #PBS -q debug-c
    #PBS -l select=1
    #PBS -l walltime=30:00
    #PBS -W group_list=group1
    #PBS -j oe
    #------- Program execution -------#
    cd ${PBS_O_WORKDIR}
    ./a.out
    ```
*   **単一ノードスレッド並列ジョブ (OpenMP)**: (例: 64スレッド)
    ```bash
    #!/bin/sh
    #------ qsub option --------#
    #PBS -q debug-c
    #PBS -l select=1:ompthreads=64
    #PBS -l walltime=01:00:00
    #PBS -W group_list=group1
    #PBS -j oe
    #------- Program execution -------#
    cd ${PBS_O_WORKDIR}
    ./a.out
    ```
    OpenMP並列化されたプログラムのスレッド数は、`OMP_NUM_THREADS` 環境変数で指定できます。
*   **単一ノードMPI並列ジョブ**: (例: ノード内112プロセス)
    ```bash
    #!/bin/sh
    #------ qsub option --------#
    #PBS -q debug-c
    #PBS -l select=1:mpiprocs=112
    #PBS -l walltime=15:00
    #PBS -W group_list=group1
    #PBS -j oe
    #------- Program execution -------#
    cd ${PBS_O_WORKDIR}
    mpiexec.hydra ./a.out
    ```
    `-l select` オプションの `mpiprocs` には、ノード内のMPI並列数を指定します。Miyabi-Cでは `mpiexec.hydra` コマンドを使用します。
*   **複数ノードMPI並列ジョブ**: (例: 2ノード、ノード内48プロセス)
    ```bash
    #!/bin/sh
    #------ qsub option --------#
    #PBS -q regular-c
    #PBS -l select=2:mpiprocs=48
    #PBS -l walltime=12:00:00
    #PBS -W group_list=group1
    #PBS -j oe
    #------- Program execution -------#
    cd ${PBS_O_WORKDIR}
    mpiexec.hydra ./a.out
    ```
    `-l select` オプションの `mpiprocs` にはノード内MPI並列数を指定し、**総MPI並列数 = mpiprocs指定数 x 要求ノード数** となるように指定してください。
*   **単一/複数ノードハイブリッド並列ジョブ (MPI+OpenMP)**: (例: 単一ノード、4プロセス、1プロセスあたり8スレッド)
    ```bash
    #!/bin/sh
    #------ qsub option --------#
    #PBS -q short-c
    #PBS -l select=1:mpiprocs=4:ompthreads=8
    #PBS -l walltime=04:00:00
    #PBS -W group_list=group1
    #PBS -j oe
    #------- Program execution -------#
    cd ${PBS_O_WORKDIR}
    mpiexec.hydra ./a.out
    ```
    ハイブリッド並列ジョブの場合、`-l select` オプションで `mpiprocs` と `ompthreads` の両方の指定が必要です。

### 5. ジョブの投入と実行状態の確認

作成したジョブスクリプトファイル（例: `jobscript.sh`）を `qsub` コマンドで投入します。

例: `qsub jobscript.sh`

ジョブが投入されると、ジョブ管理システムによって計算資源の空き状況などを元に自動で実行が開始されます。

投入したジョブの状態は `qstat` コマンドで確認できます。

例: `qstat`

`qstat` コマンドは、現在実行中または実行待ちのジョブ状態を表示します。キューの状態、ノード/MIGインスタンスの使用状況 (`--rscuse`, `--nodeuse`, `--miguse` オプション使用時)、およびジョブ本数/ノード制限 (`--limit` オプション使用時) も表示可能です。

インタラクティブジョブを実行する場合は、`qsub -I` コマンドにジョブ投入オプションを付与して実行します。Miyabi-Cのインタラクティブジョブでは、自動でIntel OneAPIがロードされます。

例 (Miyabi-C インタラクティブジョブ逐次実行):
```bash
[username@miyabi-c1 work]$ qsub -I -l select=1 -W group_list=group1 \
-q interact-c -l walltime=01:00:00
qsub: waiting for job 123456.opbs to start
qsub: job 123456.opbs ready
[username@mc001 ~]$ cd ${PBS_O_WORKDIR}
[username@mc001 work]$ icc hello_world.c
[username@mc001 work]$ ./a.out
Hello world
[username@mc001 work]$ exit
logout
qsub: job 123456.opbs completed
```

### 6. その他

ブラウザベースでの利用環境として、Open OnDemandが提供されています。利用支援ポータルからアクセスでき、仮想デスクトップ環境やJupyterLabなどが利用可能です。Open OnDemandから投入したジョブも、アクティブなジョブ確認機能で確認および削除が可能です。

コンテナでの実行環境として、singularityを利用することも可能です。

ご不明な点がございましたら、利用支援ポータル等で提供されているFAQを参照するか、問い合わせ先に連絡してください。最新の利用手引書は利用支援ポータルより入手可能です。