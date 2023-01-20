configfile: "config.yaml"

run = config.get("run", {})
RDIR = run["name"] + "/" if run.get("name") else ""

rule build_and_validate_all:
    input:
        "validation/"+RDIR+"installed_generation_capacity.xlsx",    # validate_database
        "validation/"+RDIR+"psse_cross_border_flow.xlsx",           # validate_psse

rule build_psse_network:
    name: "Build PSS/E"
    input:
        database="database/" + RDIR + "nordics.sqlite",
    output:
        loc="psse/"+RDIR+"nordics.loc",
        sav="psse/"+RDIR+"nordics.sav",
        sld="psse/"+RDIR+"nordics.sld"
    log:
        "logs/"+RDIR+"build_psse_nordics.log",
    script:
        "scripts/build_psse_network.py"

rule validate_database:
    name: "Validate Database"
    input:
        database="database/" + RDIR + "nordics.sqlite",
        entsoe_capacities="data/" + RDIR + "entsoe-transparency/installed_capacity.csv",
        entsoe_generation="data/" + RDIR + "entsoe-transparency/actual_generation.csv"
    output:
        installed_generation_capacity="validation/"+RDIR+"installed_generation_capacity.xlsx",
        actual_generation="validation/"+RDIR+"actual_generation.xlsx"
    log:
        "logs/" + RDIR + "validate_database.log",
    script:
        "scripts/validate_database.py"

rule validate_psse:
    name: "Validate PSS/E"
    input:
        sav="psse/"+RDIR+"nordics.sav",
        cross_border_flows="data/" + RDIR + "entsoe-transparency/cross_border_flow.csv",
    output:
        cross_border_flow="validation/"+RDIR+"psse_cross_border_flow.xlsx"
    log:
        "logs/" + RDIR + "validate_psse.log",
    script:
        "scripts/validate_psse.py"

rule retrieve_installed_capacity:
    name: "Retrieve installed capacity"
    output:
        installed_capacity="data/" + RDIR + "entsoe-transparency/installed_capacity.csv"
    log:
        "logs/retrieve_installed_capacity.log",
    script:
        "scripts/retrieve_installed_capacity.py"

rule retrieve_actual_generation:
    name: "Retrieve actual generation"
    output:
        actual_generation="data/" + RDIR + "entsoe-transparency/actual_generation.csv"
    log:
        "logs/" + RDIR + "retrieve_actual_generation.log"
    script:
        "scripts/retrieve_actual_generation.py"

rule retrieve_load:
    name: "Retrieve load"
    output:
        load="data/" + RDIR + "entsoe-transparency/load.csv"
    log:
        "logs/" + RDIR + "retrieve_load.log"
    script:
        "scripts/retrieve_load.py"

rule retrieve_cross_border_flow:
    output:
        cross_border_flows="data/" + RDIR + "entsoe-transparency/cross_border_flow.csv",
    log:
        "logs/retrieve_cross_border_flow.log"
    script:
        "scripts/retrieve_cross_border_flow.py"

# cross-border flows and actual generation
rule entsoe_to_sqlite:
        name: "Build database (ENTSO-E)"
        input:
            actual_generation="data/" + RDIR + "entsoe-transparency/actual_generation.csv",
            cross_border_flows="data/" + RDIR + "entsoe-transparency/cross_border_flow.csv",
            database="database/raw/" + RDIR + "nordics_without_entsoe.sqlite",
        output:
            database = "database/" + RDIR + "nordics.sqlite",
        log:
            "logs/" + RDIR + "entsoe_to_sqlite.log",
        script:
            "scripts/entsoe_to_sqlite.py"

rule pypsa_to_sqlite:
    name: "Build database (PyPSA)"
    input:
        network="input/" + RDIR + "elec.nc",
        load="data/" + RDIR + "entsoe-transparency/load.csv",
    output:
        database="database/raw/" + RDIR + "nordics_without_entsoe.sqlite",
        database_raw="database/raw/" + RDIR + "nordics.sqlite",
    log:
        "logs/" + RDIR + "pypsa_to_sqlite.log",
    script:
        "scripts/pypsa_to_sqlite.py"