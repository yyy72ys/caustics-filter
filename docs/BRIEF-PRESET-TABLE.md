# BRIEF-PRESET-TABLE — 環境プリセット表モジュールの新設＋プリセット5の登録

作業ツリー: `D:\yohey2026\dsh\muse\caustics-filter\offline`
成果物: 新規モジュール `presets_env.py` ＋自己テスト `selftest_preset_env.py` ＋報告 `out\preset-table\REPORT-PRESET-TABLE.md`

## 背景

- プリセット5「深い澄んだ池・風」の組み込み（BRIEF-PRESET5）は「コード変更ゼロ」で終わった。理由：ブリーフの形（U／F／z／媒体／位相 n・μ／η／スペクトルを持つ番号付き環境プリセット表）がツリーに存在しない。唯一の環境プリセット一覧 `web/recipe.js` の `presets()` は表示用で、F・媒体ID・FF の n/μ・η・スペクトルを持てず、代用は禁止された。
- 物理値そのものは参照解決で全項目一致済み（下記の検証値）。
- 本人決定（2026-11-15）：プリセット表モジュールを作る（発注する）。

## 仕様

新規ファイル `offline/presets_env.py` を作る。

エントリの形（schema）:
```
{id: int, name: str, U: float, fetch: float, z: float,
 medium: str, phase_key: str, eta: float, spectrum: str, arm: str}
```

登録は id=5 のみ（他の id は作らない、推測で埋めない）:
```
{5, "深い澄んだ池・風", U=6.0, fetch=100.0, z=4.0,
 medium="clear", phase_key="clear", eta=1.0,
 spectrum="elfouhaily-form", arm="D"}
```

解決関数 `resolve_preset_env(id)` を置く。BRIEF-PRESET5 の参照解決と同じ経路（`run_v15_spectrum.assemble` ＋ `v12_water.preset_spectra('clear')` ＋ `run_v16_water.FF_PHASE['clear']` ＋ `recipe_v9.parse_recipe_v9`）で完全パラメータ dict を返し、少なくとも次を含む：U, fetch, z, medium_preset, medium_label, medium_load, phase_n, phase_mu, eta, spectrum, arm, omega_c, warn_extrapolated_wave_age。実装者はツリーで実際の関数名・経路を確認し、検証値と一致する経路を使うこと（経路を変えて値だけ合わせるのは禁止）。

## 自己テスト（検証値、BRIEF-PRESET5 の実測と一致させる）

`resolve_preset_env(5)` の結果が次と一致することを assert し、11項目の結果を報告に出す：

1. U = 6.0
2. fetch = 100.0
3. z = 4.0
4. medium_preset = "clear"
5. phase_n = 1.098
6. phase_mu = 3.37
7. eta = 1.0
8. spectrum = "elfouhaily-form"
9. arm = "D"
10. omega_c = 7.296464998547599（許容差 1e-6）
11. warn_extrapolated_wave_age = True

## 制約（絶対）

- 既存ファイルは一切変更しない。新規ファイルは `presets_env.py`、`selftest_preset_env.py`、`out\preset-table\` 以下の報告のみ。
- `web/recipe.js`（sha256 8444399526442601923befb46c3e495aaba86b9b6725d1c5a7ca9af4f24bebc3）と `blender/recipe-*.json` には触らない。作業前後で主要ファイルの sha256 を比較し、無変更を報告で示す。
- 値の丸め・代用・推測での穴埋めは禁止。一致しない項目があれば、合わせに行かず「一致しない」と報告する。

## 報告に含めるもの

- モジュールの置き場と schema、登録したエントリの全文
- 自己テスト 11項目の結果（PASS/FAIL と実測値）
- 既存ファイル無変更の確認（sha256 前後比較）
- 解決に使った実際の経路（モジュール名・関数名）
