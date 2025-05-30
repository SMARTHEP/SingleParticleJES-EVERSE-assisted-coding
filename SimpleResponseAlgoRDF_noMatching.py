import ROOT
import sys
import os
ROOT.xAOD.Init(); 
from xAODDataSource import Helpers

ROOT.xAOD.CaloClusterAuxContainer_v2()
ROOT.xAOD.TruthParticleAuxContainer_v1()
# Generated using Experimental Gemini 2.5 on 28/04/2025

# --- Configuration ---
# Use an environment variable or hardcode the path to your AOD file
# Make sure this file exists and is accessible!
# You can use wildcards like "path/to/my/data*.root"
# or provide a list of filenames: ["file1.root", "file2.root"]
#input_files = os.getenv("AOD_FILE_PATH", "path/to/your/AOD.pool.root")
input_files="user.cyoung.43212508.EXT0._000438.pool.root"
# The standard TTree name in ATLAS xAOD files
tree_name = "CollectionTree"
# Output file for histograms
output_filename = "analysis_output.root"

output_file = ROOT.TFile.Open(output_filename, "RECREATE")
if not output_file or output_file.IsZombie():
    print(f"ERROR: Could not open output file {output_filename}")
    sys.exit(1)

output_file.cd()
# This would be useful if we were using the full EDM, but we're just reading the xAOD
# Check if ROOT can find the xAOD classes (a basic check for asetup)
#try:
#ROOT.xAOD.CaloCalTopoClustersAux
#except AttributeError:
#     print("ERROR: Could not find ATLAS xAOD classes (e.g., xAOD.JetContainer_v1).")
#     print("       Did you forget to run 'asetup <YourAtlasRelease>' first?")
#     sys.exit(1)


# --- RDataFrame Analysis ---

# Enable Implicit Multi-Threading (recommended for performance)
# ROOT will automatically use multiple cores if available.
ROOT.EnableImplicitMT()
print(f"Implicit Multi-Threading enabled using {ROOT.GetThreadPoolSize()} threads.")

# Create the RDataFrame
# It points to the TTree 'CollectionTree' in the specified input file(s)
print(f"Processing file(s): {input_files}")
try:
    df = Helpers.MakexAODDataFrame(input_files) #ROOT.RDataFrame(tree_name, input_files)
except Exception as e:
    print(f"ERROR: Failed to create RDataFrame. Check file path and integrity.")
    print(e)
    sys.exit(1)

# Get the initial number of events (triggers the first event loop pass if no other action taken yet)
n_total_events = df.Count().GetValue()
if n_total_events == 0:
    print("ERROR: Input file(s) contain no events in the TTree '{tree_name}'.")
    sys.exit(1)
print(f"Total events in TTree: {n_total_events}")

#for debugging purposes
#print("Getting column names:")
#print(df.GetColumnNames())


# --- Define Variables and Apply Filters ---

cluster_collection = "CaloCalTopoClusters"
particle_collection = "TruthParticles"

# Define a column 'nClusters' holding the size of the cluster collection
df_clusters_and_particles = df.Define("nClusters", f"{cluster_collection}.size()")
df_clusters_and_particles = df_clusters_and_particles.Define("nParticles", f"{particle_collection}.size()")

# Filter events: require at least 1 particle and 1 cluster in the event
df_filtered = df_clusters_and_particles.Filter("nClusters >= 1", "At least 1 cluster")
df_filtered = df_filtered.Filter("nParticles >= 1", "At least 1 particle")

# Filter more: make sure that there is at least one particle that we want (a pion, PDGID=211)
df_filtered = df_filtered.Define("pdgId_lead",f"{particle_collection}.at(0)->pdgId()")
df_filtered_PDGID = df_filtered.Filter("pdgId_lead==211","At least 1 pion")
# Define the pT of the leading cluster (assuming clusters are pt-sorted, common but check!)
# ATLAS stores energies/momenta in MeV, so divide by 1000 for GeV.
# Accessing elements of the xAOD container directly:
df_defined = df_filtered_PDGID.Define("leading_cluster_e", f"{cluster_collection}.at(0)->rawE() / 1000.0")


# Define code for response, without including matching
# https://lup.lub.lu.se/luur/download?func=downloadFile&recordOId=8996766&fileOId=8996773, page 32:
# This algorithm iterates over all the events in the NTuple, and for each event iterates over all the truth particles. 
# For each truth particle it geometrically matches all reconstructed clusters within DeltaR = 0.2  to that truth particle. 
# The four-vectors of the clusters matched are then summed into one four-vector. 
# Finally, the pT response is then calculated from the pT of the combined cluster four-vector and the pT of the truth particle. 
# What we are doing here: taking the leading pion as reference and the leading cluster as probe
# This is justified because there is a reasonable correlation between the energy, eta, phi of these

response_code = f"""

    // Check if there are indeed at least 2 clusters (filter should ensure this, but defensive coding)
    if ({cluster_collection}.size() < 2) {{
        return -1.0f; // Return a dummy value
    }}
    
    auto& clusters = {cluster_collection}; // Reference for convenience
    auto& particles = {particle_collection}; // 
    
    // Create Lorentz vectors for the leading cluster and particle (assuming pt-sorted)
    // Since LorentzVector wants pT and not energy, we'll give it pT

    float cluster_pt = sqrt(clusters.at(0)->rawE()*clusters.at(0)->rawE() - clusters.at(0)->rawM()*clusters.at(0)->rawM())/TMath::CosH(clusters.at(0)->rawEta());

    ROOT::Math::PtEtaPhiMVector cluster_lv(cluster_pt, clusters.at(0)->rawEta(), clusters.at(0)->rawPhi(), clusters.at(0)->rawM());
    ROOT::Math::PtEtaPhiMVector particle_lv(particles.at(0)->pt(), particles.at(0)->eta(), clusters.at(0)->phi(), clusters.at(0)->m());

    
    // pT response 
    float response = cluster_lv.Pt()/particle_lv.Pt();
    return response;

"""

df_defined = df_defined.Define("response", response_code)
df_defined = df_defined.Define("lead_particle_pt", f"{particle_collection}.at(0)->pt()/1000.")
df_defined = df_defined.Define("lead_particle_eta", f"{particle_collection}.at(0)->eta()")

# --- Book Histograms ---
# Histograms are booked here, but only filled when the event loop is triggered later.

# Histogram model: (name, title; x-axis label; y-axis label, nbins, xlow, xhigh)
h_nClusters_model = ROOT.RDF.TH1DModel("h_nClusters", "Number of Clusters;N_{Clusters};Events", 20, -0.5, 19.5)
h_nParticles_model = ROOT.RDF.TH1DModel("h_nParticles", "Number of Particles;N_{Particles};Events", 20, -0.5, 19.5)
h_leading_cluster_e_model = ROOT.RDF.TH1DModel("h_leading_cluster_e", "Leading Cluster E;E^{lead cluster} [GeV];Events", 100, 0, 500)
h_PDGIDs_model = ROOT.RDF.TH1DModel("h_PDGIDs", "Leading Cluster PDGID;PDGID^{lead cluster};Events", 1000, 0, 1000)
h_inclusive_response_model = ROOT.RDF.TH1DModel("h_inclusive_response", "Response (inclusive); p_{T,cluster}/p_{T,particle};Events", 100, 0, 3)
h_3d_response_model = ROOT.RDF.TH3DModel("h_response", "Response (binned); p_{T,cluster}/p_{T,particle},p_{T,particle},eta;Events", 50,0,2, 100,0,500, 60,-3,3)

# Book the histograms using the defined columns
h_nClusters = df_clusters_and_particles.Histo1D(h_nClusters_model, "nClusters") # Booked on df_clusters_and_particles (before filter)
h_nParticles = df_clusters_and_particles.Histo1D(h_nParticles_model, "nParticles") # Booked on df_clusters_and_particles (before filter)
h_PDGIDs = df_filtered_PDGID.Histo1D(h_PDGIDs_model, "pdgId_lead") # Booked on df_clusters_and_particles (before filter)
h_leading_cluster_e = df_defined.Histo1D(h_leading_cluster_e_model, "leading_cluster_e") # Booked on df_defined (after filter)
h_inclusive_response = df_defined.Histo1D(h_inclusive_response_model, "response") # Booked on df_defined (after filter)
h_3d_response = df_defined.Histo3D(h_3d_response_model, "response","lead_particle_pt","lead_particle_eta") # Booked on df_defined (after filter)

h_nClusters.SetDirectory(output_file);
h_nParticles.SetDirectory(output_file);
h_PDGIDs.SetDirectory(output_file);
h_leading_cluster_e.SetDirectory(output_file);
h_inclusive_response.SetDirectory(output_file);
h_3d_response.SetDirectory(output_file);

# --- Trigger Execution and Save Output ---

# Accessing the histograms' values or using Snapshot triggers the actual event loop.
# RDataFrame processes events lazily only when results are requested.

print("Event loop running...")

# You can also generate a report on filter efficiencies, etc.
report = df_defined.Report()

# Optionally, save selected data to a new, smaller ROOT file (a "Snapshot")
# Specify the columns you want to keep.
# Note: Snapshotting complex xAOD objects directly might require extra steps or helper libraries.
#       It's often easier to snapshot derived primitive types.
# columns_to_save = ROOT.std.vector['string'](["runNumber", "eventNumber", "leading_cluster_pt", "response"])
# df_defined.Snapshot("SelectedEventsTree", "snapshot_output.root", columns_to_save)
# print("Snapshot file created: snapshot_output.root")


# Save histograms to the output file
print(f"Saving histograms to: {output_filename}")

output_file.cd()

# Calling GetValue() here ensures the loop ran and histograms are filled *before* writing
h_nClusters_ptr = h_nClusters.GetValue()
h_nParticles_ptr = h_nParticles.GetValue()
h_PDGIDs_ptr = h_PDGIDs.GetValue()
h_leading_cluster_e_ptr = h_leading_cluster_e.GetValue()
h_inclusive_response_ptr = h_inclusive_response.GetValue()
h_3d_response_ptr = h_3d_response.GetValue()

h_nClusters_ptr.Write()
h_nParticles_ptr.Write()
h_PDGIDs_ptr.Write()
h_leading_cluster_e_ptr.Write()
h_inclusive_response_ptr.Write()
h_3d_response_ptr.Write()

output_file.Close()
print("Histograms saved.")

# Print the report
print("\n--- Analysis Report ---")
report.Print()
print("---------------------\n")

print("Analysis finished successfully.")
