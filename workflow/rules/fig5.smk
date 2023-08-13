def script_path(x):
    return f"workflow/scripts/{x}"


rule all:
    input:
        "figures/figure5/simulation_ins_mouse.png",
        "figures/figure5/conditions_mouse.png",


rule simulation:
    output:
        mouse="figures/figure5/simulation_ins_mouse.png",
        human="figures/figure5/simulation_ins_human.png",
        turtle="figures/figure5/simulation_ins_turtle.png",
    script:
        script_path("fig5_simulation.py")


rule conditions:
    output:
        mouse="figures/figure5/conditions_mouse.png",
        human="figures/figure5/conditions_human.png",
        turtle="figures/figure5/conditions_turtle.png",
    script:
        script_path("fig5_conditions.py")
