# L24 next steps

L24 was intentionally inactive for this release.

1. Confirm/resume the paused L24/110_A2 and L24/111_A2 determinant scans.
2. Verify the cache metadata Ecm window and corrected raw-reader provenance.
3. Generate determinant-family plots and candidate brackets for each L24 irrep.
4. Manually label every candidate; only explicit `user_label=true` rows may be promoted.
5. Resolve L24/100_A2 bracket 36 before including that sector.
6. Build `accepted_windows_L24_user_validated_*` and expanded projected-basis hot-window caches.
7. Run row-coverage and metadata/hash audits.
8. Run one-FCN tests with `fallback_full_scan=0` before any minimization.
9. Run combined L20+L24 fits only after sector-level classifier and model-root counts pass.

Do not run cachegen, use false/blank windows, or enable original-cache fallback.
