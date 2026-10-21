# Large-file and cache policy

The five final L20 hot-window projected-basis binaries are included in the
full working zip but ignored by normal Git because three individual files are
over 100 MB. The original full raw F3inv/Vsel cache root is never copied.

To restore the hot-window payloads from the zip, extract the archive over the
release tree and verify:

```bash
find output/hot_windows_L20_user_validated_0p36 -name '*.bin' -printf '%s %p\n'
```

The JSON metadata, validation CSVs, accepted-window hashes, reports, and all
fit outputs are tracked in Git. A fresh fit still requires the external raw
cache root documented in the input configurations.
