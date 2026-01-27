def export_parameter_tables():
    # Occurrence table
    occ = raw_conc_df.copy()
    occ['Distribution'] = 'Lognormal'
    occ['Sigma_assumed'] = 0.5  # if you override per-pathogen, join those here
    occ = occ.rename(columns={
        'Pathogen': 'Pathogen',
        'Matrix': 'Matrix',
        'Raw_Concentration': 'Geometric_mean_conc'  # interpret as GM; or rename accordingly
    })
    occ[['Pathogen','Matrix','Distribution','Geometric_mean_conc','Sigma_assumed']].to_csv(
        'occurrence_parameters_table.csv', index=False
    )

    # Treatment table (flatten piece-wise)
    rows = []
    for proc, d in TREATMENT_EFFECTIVENESS_DIST.items():
        for ptype, profiles in d.items():
            for prof in profiles:
                rows.append({
                    'Process': proc,
                    'Pathogen_type': ptype,
                    'Distribution': 'Truncated normal',
                    'LR_mean': prof['mean'],
                    'LR_sd': prof['sd'],
                    'Max_conc_threshold': prof['max_conc']
                })
    pd.DataFrame(rows).to_csv('treatment_parameters_table.csv', index=False)
    print("Exported: occurrence_parameters_table.csv, treatment_parameters_table.csv")
