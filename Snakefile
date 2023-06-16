configfile: "config.yaml"

run = config.get("run", {})
RDIR = run["name"] + "/" if run.get("name") else ""
RDIR_INPUT = run["input"] + "/" if run.get("input") else ""

rule build_and_validate_all:
    input:
        "database/"+RDIR+"nordics.sqlite",                                      # build_database
        "psse/"+RDIR+"nordics.sav",                                             # build_psse
        "validation/"+RDIR+"installed_generation_capacity_per_country.xlsx",
        "validation/"+RDIR+"balance.xlsx",                                      # validate_database
        "psse/"+RDIR+"nordics.sld",                                             # build_psse_slider
        "validation/"+RDIR+"psse_cross_border_flow.xlsx",                       # validate_psse

rule build_psse_network:
    name: "Build PSS/E"
    input:
        database="database/" + RDIR + "nordics.sqlite",
        remove_lines="grid-change/remove_lines.xlsx",
        #inactive_lines="manual-data/" + RDIR + "inactive-lines.xlsx",
        remove_buses = "grid-change/remove-bus.xlsx",
        remove_transformers = "grid-change/remove_transformers.xlsx",
    output:
        sav="psse/"+RDIR+"nordics.sav",
    log:
        "logs/"+RDIR+"build_psse_nordics.log",
    script:
        "scripts/build_psse_network.py"


rule build_psse_slider:
    name: "Build PSS/E slider"
    input:
        database="database/" + RDIR + "nordics.sqlite",
        sav="psse/" + RDIR + "nordics.sav",
    output:
        loc="psse/"+RDIR+"nordics.loc",
        sld="psse/"+RDIR+"nordics.sld"
    log:
        "logs/"+RDIR+"build_psse_slider.log",
    script:
        "scripts/build_psse_slider.py"

rule validate_database:
    name: "Validate Database"
    input:
        database="database/" + RDIR + "nordics_raw.sqlite",
        entsoe_cross_border_flow="data/" + RDIR + "entsoe-transparency/cross_border_flow.csv",
        entsoe_capacities="data/" + RDIR + "entsoe-transparency/installed_capacity.csv",
        entsoe_generation="data/" + RDIR + "entsoe-transparency/actual_generation.csv",
        entsoe_load="data/" + RDIR + "entsoe-transparency/load.csv",
    output:
        generation_capacity_country="validation/"+RDIR+"installed_generation_capacity_per_country.xlsx",
        generation_capacity_bidding_zone="validation/"+RDIR+"installed_generation_capacity_per_bidding_zone.xlsx",
        actual_generation_country="validation/"+RDIR+"actual_generation_per_country.xlsx",
        actual_generation_bidding_zone="validation/"+RDIR+"actual_generation_per_bidding_zone.xlsx",
        balance="validation/"+RDIR+"balance.xlsx",
    log:
        "logs/" + RDIR + "validate_database.log",
    script:
        "scripts/validate_database.py"

rule merge_databases:
    input:
        bus_info_db_path="database/" + RDIR + "nordics.sqlite",
        bus_location_db_path="example.db",
    output:
        "database/" + RDIR + "merged_db.sqlite"
    shell:
        "python merge_databases.py {input.bus_info_db_path} {input.bus_location_db_path} {output.merged_db_path}"


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
    name: "Retrieve cross-border flow"
    output:
        cross_border_flows="data/" + RDIR + "entsoe-transparency/cross_border_flow.csv",
    log:
        "logs/retrieve_cross_border_flow.log"
    script:
        "scripts/retrieve_cross_border_flow.py"

# cross-border flows and actual generation
rule build_database:
        name: "Build database"
        input:
            network="input/" + RDIR_INPUT + "elec.nc",
            actual_generation="data/" + RDIR + "entsoe-transparency/actual_generation.xlsx",
            cross_border_flows="data/" + RDIR + "entsoe-transparency/cross_border_flow.csv",
            load = "data/" + RDIR + "entsoe-transparency/load.csv",
            bus_names = "grid-change/bus-names.xlsx",
            manual_loads = "manual-data/" + RDIR + "loads.xlsx",
            manual_shunts = "grid-change/shunts.xlsx",
            manual_generators = "manual-data/" + RDIR + "generators.xlsx",
            manual_buses = "grid-change/buses.xlsx",
            manual_lines = "grid-change/lines.xlsx",
            generators_case = "manual-data/" + RDIR + "generators_case.xlsx",
            case_lines = "manual-data/" + RDIR + "case_lines.xlsx",
        output:
            database = "database/" + RDIR + "nordics.sqlite",
            database_raw = "database/" + RDIR + "nordics_raw.sqlite",
        log:
            "logs/" + RDIR + "build_database.log",
        script:
            "scripts/build_database.py"

rule load_duration_curve:
        name: "Load duration curve"
        input:
            network = "input/" + RDIR_INPUT + "elec.nc",
        output:
            load_duration_curve="data/" + RDIR + "network_load_duration_curve.png",
        log:
            "logs/" + RDIR + "load_duration_curve.log"
        script:
            "scripts/load_duration_curve.py"

# rule merge_buses_rule:
#        input:
#            nordics_db="database/" + RDIR + "nordics.sqlite",
#            example_db="example.sqlite"
#        output:
#            merged_db="database/" + RDIR + "nordics.sqlite"
#        run:
#            merge_buses(nordics_db, example_db, output.merged_db)
