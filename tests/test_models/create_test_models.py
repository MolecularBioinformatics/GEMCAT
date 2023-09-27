import cobra

# Mini

a = cobra.Metabolite(id="A", name="A")
b = cobra.Metabolite(id="B", name="B")
c = cobra.Metabolite(id="C", name="C")
d = cobra.Metabolite(id="D", name="D")
mini_mets = [a, b, c, d]

r1 = cobra.Reaction(id="R1", name="R1", lower_bound=0, upper_bound=1)
r2 = cobra.Reaction(id="R2", name="R2", lower_bound=0, upper_bound=1)
r3 = cobra.Reaction(id="R3", name="R3", lower_bound=0, upper_bound=1)
r4 = cobra.Reaction(id="R4", name="R4", lower_bound=0, upper_bound=1)
mini_rxns = [r1, r2, r3, r4]

for rxn in mini_rxns:
    rxn.compartment = "c"
for met in mini_mets:
    met.compartment = "c"

r1.add_metabolites({a: -1, b: 1})
r2.add_metabolites({a: -1, c: 1})
r3.add_metabolites({b: -1, d: 1})
r4.add_metabolites({c: -1, d: 1})

g1 = cobra.Gene(id="G1")
g2 = cobra.Gene(id="G2")
g3 = cobra.Gene(id="G3")
g4 = cobra.Gene(id="G4")

r1.gene_reaction_rule = "G1"
r2.gene_reaction_rule = "G2"
r3.gene_reaction_rule = "G3"
r4.gene_reaction_rule = "G4"

mini = cobra.Model("test_mini", "test_mini")
mini.add_reactions(mini_rxns)

cobra.io.write_sbml_model(mini, "mini.xml")
cobra.io.save_json_model(mini, "mini.json")

# Mini reversible

mini.reactions.get_by_id("R1").lower_bound = -1
cobra.io.write_sbml_model(mini, "mini_reversible.xml")
cobra.io.save_json_model(mini, "mini_reversible.json")

# # Mini redox

mini.reactions.get_by_id("R1").lower_bound = 0

e = cobra.Metabolite("E", compartment="c")
f = cobra.Metabolite("F", compartment="c")
r3 = mini.reactions.get_by_id("R3")
r4 = mini.reactions.get_by_id("R4")
r3.add_metabolites({e: -1.0, f: 1.0})
r4.add_metabolites({e: 1.0, f: -1.0})
r4.metabolites

cobra.io.write_sbml_model(mini, "mini_redox.xml")
cobra.io.save_json_model(mini, "mini_redox.json")

# # Mini complex GPR

g1 = cobra.Gene(id="G1")
g2 = cobra.Gene(id="G2")
g3 = cobra.Gene(id="G3")
g4 = cobra.Gene(id="G4")
g5 = cobra.Gene(id="G5")
g6 = cobra.Gene(id="G6")
g7 = cobra.Gene(id="G7")
g8 = cobra.Gene(id="G8")

mini.reactions.R3.gene_reaction_rule = "G3 or (G4 and G5 and G6) or (G7 and G8)"

cobra.io.write_sbml_model(mini, "mini_complex_gpr.xml")
cobra.io.save_json_model(mini, "mini_complex_gpr.json")

g1 = cobra.Gene(id="0001.1")
g2 = cobra.Gene(id="0002.2")
g3 = cobra.Gene(id="0003.3")
g4 = cobra.Gene(id="0004.4")
g5 = cobra.Gene(id="0005.5")
g6 = cobra.Gene(id="0006.6")
g7 = cobra.Gene(id="0007.7")
g8 = cobra.Gene(id="0008.8")

mini.reactions.R1.gene_reaction_rule = "0001.1"
mini.reactions.R2.gene_reaction_rule = "0002.2"
mini.reactions.R3.gene_reaction_rule = (
    "0003.3 or (0004.4 and 0005.5 and 0006.6) or (0007.7 and 0008.8)"
)
mini.reactions.R4.gene_reaction_rule = "0004.4"

cobra.io.write_sbml_model(mini, "mini_complex_gpr_gid_version.xml")
cobra.io.save_json_model(mini, "mini_complex_gpr_gid_version.json")

# # Highly connected

hc = cobra.Model("Highly connected")

a = cobra.Metabolite("A", compartment="c")
b = cobra.Metabolite("B", compartment="c")
c = cobra.Metabolite("C", compartment="c")
d = cobra.Metabolite("D", compartment="c")
e = cobra.Metabolite("E", compartment="c")
f = cobra.Metabolite("F", compartment="c")
g = cobra.Metabolite("G", compartment="c")

r1 = cobra.Reaction("R1")
r2 = cobra.Reaction("R2")
r3 = cobra.Reaction("R3")
r4 = cobra.Reaction("R4")
r5 = cobra.Reaction("R5")
r6 = cobra.Reaction("R6")
r7 = cobra.Reaction("R7")
r8 = cobra.Reaction("R8")
r9 = cobra.Reaction("R9")
r10 = cobra.Reaction("R10")
r11 = cobra.Reaction("R11")

r1.gene_reaction_rule = "G1"
r2.gene_reaction_rule = "G2"
r3.gene_reaction_rule = "G3"
r4.gene_reaction_rule = "G4"
r5.gene_reaction_rule = "G5"
r6.gene_reaction_rule = "G6"
r7.gene_reaction_rule = "G7"
r8.gene_reaction_rule = "G8"
r9.gene_reaction_rule = "G9"
r10.gene_reaction_rule = "G10"
r11.gene_reaction_rule = "G11"

r1.add_metabolites({a: -1, b: 1})
r2.add_metabolites({a: -1, d: 1})
r3.add_metabolites({b: -1, c: 1})
r4.add_metabolites({d: -1, e: 1})
r5.add_metabolites({b: -1, f: 1})
r6.add_metabolites({d: -1, f: 1})
r7.add_metabolites({f: -1, c: 1})
r8.add_metabolites({f: -1, e: 1})
r9.add_metabolites({c: -1, g: 1})
r10.add_metabolites({e: -1, g: 1})
r11.add_metabolites({f: -1, g: 1})

hc.add_reactions([r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11])

cobra.io.write_sbml_model(hc, "highly_connected.xml")
cobra.io.save_json_model(hc, "highly_connected.json")
