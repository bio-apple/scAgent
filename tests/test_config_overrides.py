from scagent.config import apply_performance_overrides, dask_params, gpu_params, load_config


def test_apply_performance_overrides_dask():
    cfg = apply_performance_overrides(load_config(reload=True), dask=True)
    assert dask_params(cfg)["enabled"] is True


def test_apply_performance_overrides_rapids_enables_gpu():
    cfg = apply_performance_overrides(load_config(reload=True), rapids=True)
    assert gpu_params(cfg)["enabled"] is True
    assert gpu_params(cfg)["rapids"] is True
    load_config(reload=True)
