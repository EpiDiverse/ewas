class ParameterChecks {
	static void checkParams(params) {
		assert params.input, "Please specify path to input with --input parameter"
        assert params.samples, "Please specify path to samples.tsv with --samples parameter"
        assert !(params.Gmodel && !params.SNPs), "Please specify path to variants with --SNPs parameter when running --Gmodel"
        assert !(params.GxE && !params.SNPs), "Please specify path to variants with --SNPs parameter when running --GxE"
        assert !(params.noCpG && params.noCHG && params.noCHH), "Please specify at least one methylation context for analysis!"
        assert params.coverage instanceof Integer && params.coverage >= 0, "Input coverage filter must be a non-negative integer!"
        assert params.proportion instanceof Numeric && params.proportion in 0..1, "Proportion of shared regions must be a decimal in the range of 0 and 1!"
        assert params.input_FDR instanceof Float && params.input_FDR in 0..1, "Input FDR filter must be a decimal in the range of 0 and 1!"
        assert params.output_FDR instanceof Float && params.output_FDR in 0..1, "Output FDR filter must be a decimal in the range of 0 and 1!"
        assert params.Emodel_pv instanceof Float && params.Emodel_pv in 0..1, "Emodel p-value must be a decimal in the range of 0 and 1!"
        assert params.Gmodel_pv instanceof Float && params.Gmodel_pv in 0..1, "Gmodel p-value must be a decimal in the range of 0 and 1!"
        assert params.GxE_pv instanceof Float && params.GxE_pv in 0..1, "GxE p-value must be a decimal in the range of 0 and 1!"
        assert params.max_missing instanceof Float && params.max_missing in 0..1, "--max_missing parameter must be a decimal in the range of 0 and 1!"
        assert params.mac instanceof Integer && params.mac >= 0, "--mac parameter must be a non-negative integer!"
        assert params.minQ instanceof Integer && params.minQ >= 0, "--minQ parameter must be a non-negative integer!"
        assert params.take instanceof Integer && params.take >= 0, "--take parameter must be a non-negative integer!"
        assert params.fork instanceof Integer && params.fork >= 0, "--fork parameter must be a non-negative integer!"
	}
}