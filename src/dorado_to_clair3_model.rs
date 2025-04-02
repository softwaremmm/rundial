use phf::phf_map;

static MAPPING_TABLE: phf::Map<&'static str, &'static str> = phf_map! {
    "dna_r9.4.1_450bps_sup_prom" => "r941_prom_sup_g5014",
    "dna_r10.4.1_e8.2_260bps_fast@v3.5.2" => "r1041_e82_260bps_fast_g632",
    "dna_r10.4.1_e8.2_260bps_hac@v3.5.2" => "r1041_e82_260bps_hac_g632",
    "dna_r10.4.1_e8.2_260bps_sup@v3.5.2" => "r1041_e82_260bps_sup_g632",
    "dna_r10.4.1_e8.2_400bps_fast@v3.5.2" => "r1041_e82_400bps_fast_g632",
    "dna_r10.4.1_e8.2_400bps_hac@v3.5.2" => "r1041_e82_400bps_hac_g632",
    "dna_r10.4.1_e8.2_400bps_sup@v3.5.2" => "r1041_e82_400bps_sup_g615",
    "dna_r10.4.1_e8.2_260bps_fast@v4.0.0" => "r1041_e82_260bps_hac_v400",
    "dna_r10.4.1_e8.2_260bps_hac@v4.0.0" => "r1041_e82_260bps_hac_v400",
    "dna_r10.4.1_e8.2_260bps_sup@v4.0.0" => "r1041_e82_260bps_sup_v400",
    "dna_r10.4.1_e8.2_400bps_fast@v4.0.0" => "r1041_e82_400bps_hac_v400",
    "dna_r10.4.1_e8.2_400bps_hac@v4.0.0" => "r1041_e82_400bps_hac_v400",
    "dna_r10.4.1_e8.2_400bps_sup@v4.0.0" => "r1041_e82_400bps_sup_v400",
    "dna_r10.4.1_e8.2_260bps_fast@v4.1.0" => "r1041_e82_260bps_hac_v410",
    "dna_r10.4.1_e8.2_260bps_hac@v4.1.0" => "r1041_e82_260bps_hac_v410",
    "dna_r10.4.1_e8.2_260bps_sup@v4.1.0" => "r1041_e82_260bps_sup_v410",
    "dna_r10.4.1_e8.2_400bps_fast@v4.1.0" => "r1041_e82_400bps_hac_v410",
    "dna_r10.4.1_e8.2_400bps_hac@v4.1.0" => "r1041_e82_400bps_hac_v410",
    "dna_r10.4.1_e8.2_400bps_sup@v4.1.0" => "r1041_e82_400bps_sup_v410",
    "dna_r10.4.1_e8.2_400bps_fast@v4.2.0" => "r1041_e82_400bps_hac_v420",
    "dna_r10.4.1_e8.2_400bps_hac@v4.2.0" => "r1041_e82_400bps_hac_v420",
    "dna_r10.4.1_e8.2_400bps_sup@v4.2.0" => "r1041_e82_400bps_sup_v420",
    "dna_r10.4.1_e8.2_400bps_fast@v4.3.0" => "r1041_e82_400bps_hac_v430",
    "dna_r10.4.1_e8.2_400bps_hac@v4.3.0" => "r1041_e82_400bps_hac_v430",
    "dna_r10.4.1_e8.2_400bps_sup@v4.3.0" => "r1041_e82_400bps_sup_v430",
};

/// static function which simply applies mapping table
pub fn map_model_name(model_name: &str) -> Option<&'static str> {
    MAPPING_TABLE.get(model_name).copied()
}
