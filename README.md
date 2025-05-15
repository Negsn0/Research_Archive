# Research_Archive

## 概要

原子核の有限温度における超流動相転移の研究アーカイブです。Woods-Saxon ポテンシャルと seniority pairing モデルを用いて、Sn 同位体の超流動相転移を解析しています。

## ディレクトリ構成

### メインディレクトリ

- `Presentation/`: プレゼンテーション関連ファイル
  - `fig_pdf/`: 図表ファイル
  - `0228_fix.tex`: 発表スライド（修正版）
  - `0228.tex`: 発表スライド（元版）
- `Graduation_thesis/`: 卒業論文関連ファイル
- `complete_code/`: 完成版のプログラムコード
- `reserch_try/`: 研究過程での試行錯誤コード
- `output/`: 計算結果の出力ファイル
- `参考文献/`: 参考文献資料

### プログラムファイル

- `SnEnergyLevels.jl`: Sn 同位体のエネルギー準位計算プログラム

### 論文・ドキュメント

- `main.tex`, `main_fix.tex`: 本論文（原版・修正版）
- `main_before.tex`: 論文初期版
- `FTBCS.tex`: 有限温度 BCS 理論の解説
- `Fortran.tex`: Fortran コードの説明
- `Idea.tex`: 研究アイデアのメモ

### 設定ファイル

- `.gitignore`: Git の除外設定ファイル

## 研究内容

- 対象核種: \ce{^{100-132}Sn}
- 解析手法:
  - Woods-Saxon ポテンシャルによる平均場計算
  - Seniority pairing モデルによる残留相互作用の取り扱い
  - 有限温度 BCS 理論による超流動相転移の解析

## 主な結果

- 全核種で超流動相転移を確認
- 相転移温度は$kT/\Delta_0 \sim 0.55$付近に集中
- エネルギー期待値と比熱の温度依存性を評価

## 使用言語・環境

- LaTeX (beamer クラス)
- 図表作成: PDF 形式
- Julia
- Fortran

## 参考文献

詳細は発表スライド内の参考文献リストを参照してください。
