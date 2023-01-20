configfile: "config.yaml"

CASES = config.get('snapshots').keys()
COUNTRIES = config.get('countries')

run = config.get("run", {})
RDIR = run["name"] + "/" if run.get("name") else ""

rule build_and_validate_all:
    input:
        expand("validation/"+RDIR+"{case}/installed_generation_capacity.xlsx", case=CASES),    # validate_database
        expand("psse/"+RDIR+"{case}/nordics_{case}.sav", case=CASES),                          # build_psse_network
        expand("validation/"+RDIR+"{case}/psse_cross_border_flow.xlsx", case=CASES),           # validate_psse

rule build_psse_network:
    name: "Build PSS/E"
    input:
        database="database/" + RDIR + "nordics_{case}.sqlite",
    output:
        loc="psse/"+RDIR+"{case}/nordics_{case}.loc",
        sav="psse/"+RDIR+"{case}/nordics_{case}.sav",
        sld="psse/"+RDIR+"{case}/nordics_{case}.sld"
    log:
        "logs/"+RDIR+"build_psse_nordics_{case}.log",
    script:
        "scripts/build_psse_network.py"

rule validate_database:
    name: "Validate Database"
    input:
        database="database/" + RDIR + "nordics_{case}.sqlite",
        entsoe_capacities="data/" + RDIR + "entsoe-transparency/installed_capacity.csv",
        entsoe_generation="data/" + RDIR + "entsoe-transparency/actual_generation_{case}.csv"
    output:
        installed_generation_capacity="validation/"+RDIR+"{case}/installed_generation_capacity.xlsx",
        actual_generation="validation/"+RDIR+"{case}/actual_generation.xlsx"
    log:
        "logs/" + RDIR + "validate_database_{case}.log",
    script:
        "scripts/validate_database.py"

rule validate_psse:
    name: "Validate PSS/E"
    input:
        sav="psse/"+RDIR+"{case}/nordics_{case}.sav",
        cross_border_flows="data/" + RDIR + "entsoe-transparency/cross_border_flow_{case}.csv",
    output:
        cross_border_flow="validation/"+RDIR+"{case}/psse_cross_border_flow.xlsx"  # TODO: implement
    log:
        "logs/" + RDIR + "validate_psse_{case}.log",
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
    input:
        snapshots="data/" + RDIR + "snapshots.csv"
    output:
        actual_generation="data/" + RDIR + "entsoe-transparency/actual_generation_{case}.csv"
    log:
        "logs/" + RDIR + "retrieve_actual_generation_{case}.log"
    script:
        "scripts/retrieve_actual_generation.py"

rule retrieve_load:
    name: "Retrieve load"
    input:
        snapshots="data/" + RDIR + "snapshots.csv"
    output:
        load="data/" + RDIR + "entsoe-transparency/load_{case}.csv"
    log:
        "logs/" + RDIR + "retrieve_load_{case}.log"
    script:
        "scripts/retrieve_load.py"

rule retrieve_cross_border_flow:
    name: "Retrieve cross-border flow"
    input:
        snapshots = "data/" + RDIR + "snapshots.csv",
    output:
        cross_border_flows="data/" + RDIR + "entsoe-transparency/cross_border_flow_{case}.csv",
    log:
        "logs/retrieve_cross_border_flow_{case}.log"
    script:
        "scripts/retrieve_cross_border_flow.py"

# cross-border flows and actual generation
rule entsoe_to_sqlite:
        name: "Build database (ENTSO-E)"
        input:
            actual_generation="data/" + RDIR + "entsoe-transparency/actual_generation_{case}.csv",
            cross_border_flows="data/" + RDIR + "entsoe-transparency/cross_border_flow_{case}.csv",
            snapshots="data/" + RDIR + "snapshots.csv",
            database="database/raw/" + RDIR + "nordics_{case}_without_entsoe.sqlite",
        output:
            database = "database/" + RDIR + "nordics_{case}.sqlite",
        log:
            "logs/" + RDIR + "{case}_entsoe_to_sqlite.log",
        script:
            "scripts/entsoe_to_sqlite.py"


if config.get('source') == 'pypsa':
    rule pypsa_snapshots:
        name: "Snapshots"
        input:
            network = "input/" + RDIR + "elec.nc",
        output:
            load_duration_curve="data/" + RDIR + "network_load_duration_curve.png",
            snapshots="data/" + RDIR + "snapshots.csv",
        log:
            "logs/" + RDIR + "pypsa_snapshots.log"
        script:
            "scripts/pypsa_snapshots.py"

    rule pypsa_to_sqlite:
        name: "Build database (PyPSA)"
        input:
            network="input/" + RDIR + "elec.nc",
            load="data/" + RDIR + "entsoe-transparency/load_{case}.csv",
            snapshots="data/" + RDIR + "snapshots.csv",
        output:
            database="database/raw/" + RDIR + "nordics_{case}_without_entsoe.sqlite",
            database_raw="database/raw/" + RDIR + "nordics_{case}.sqlite",
        log:
            "logs/" + RDIR + "{case}_pypsa_to_sqlite.log",
        script:
            "scripts/pypsa_to_sqlite.py"