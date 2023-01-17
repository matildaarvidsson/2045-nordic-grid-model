configfile: "config.yaml"

CASES = config.get('snapshots').keys()
COUNTRIES = config.get('countries')

run = config.get("run", {})
RDIR = run["name"] + "/" if run.get("name") else ""

rule build_and_validate_all:
    input:
        expand("validation/"+RDIR+"{case}/installed_generation_capacity.xlsx", case=CASES),    # validate_database
        expand("psse/"+RDIR+"{case}/network_{case}.sav", case=CASES),                          # build_psse_network
        # expand("validation/"+RDIR+"{case}/psse_validation_report.xlsx", case=CASES),         # validate_psse


rule build_psse_network:
    name: "Build PSS/E"
    input:
        database="database/" + RDIR + "network_{case}.sqlite",
    output:
        loc="psse/"+RDIR+"{case}/network_{case}.loc",
        sav="psse/"+RDIR+"{case}/network_{case}.sav",
        sld="psse/"+RDIR+"{case}/network_{case}.sld"
    log:
        "logs/"+RDIR+"build_psse_network_{case}.log",
    script:
        "scripts/build_psse_network.py"

rule validate_database:
    name: "Validate Database"
    input:
        database="database/" + RDIR + "network_{case}.sqlite",
        entsoe_capacities="data/entsoe-transparency/installed_capacity.csv",
        entsoe_generation="data/entsoe-transparency/actual_generation_{case}.csv"
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
        sav="psse/"+RDIR+"{case}/network_{case}.sav"
    output:
        report="validation/"+RDIR+"{case}/psse_validation_report.xlsx"  # TODO: implement
    log:
        "logs/" + RDIR + "validate_psse_{case}.log",
    script:
        "scripts/validate_psse.py"

rule retrieve_installed_capacity:
    name: "Retrieve installed capacity"
    output:
        installed_capacity="data/entsoe-transparency/installed_capacity.csv"
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

rule retrieve_cross_border_flow:
    name: "Retrieve cross-border flow"
    output:
        cross_border_flows="data/entsoe-transparency/cross_border_flow.csv"
    log:
        "logs/retrieve_cross_border_flow.log"
    script:
        "scripts/retrieve_cross_border_flow.py"

# cross-border flows and actual generation
rule entsoe_to_sqlite:
        name: "Build database (ENTSO-E)"
        input:
            actual_generation="data/entsoe-transparency/actual_generation_{case}.csv",
            cross_border_flows="data/entsoe-transparency/cross_border_flow.csv",
            snapshots="data/" + RDIR + "snapshots.csv",
            database="database/raw/" + RDIR + "network_{case}_without_entsoe.sqlite",
        output:
            database = "database/" + RDIR + "network_{case}.sqlite",
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
            snapshots="data/" + RDIR + "snapshots.csv",
        output:
            database="database/raw/" + RDIR + "network_{case}_without_entsoe.sqlite",
            database_raw="database/raw/" + RDIR + "network_{case}.sqlite",
        log:
            "logs/" + RDIR + "{case}_pypsa_to_sqlite.log",
        script:
            "scripts/pypsa_to_sqlite.py"