● Great, plan locked in: Metropolis-rate KMC, BulkMonomer-only v1, time-calibrated, sim-time sampling + wall-clock cap. Here's the phased plan.
                                                                                                                                                               
  ---                                                                                                         
  Design decisions (locked before coding)                                                                                                                                                              
  - Rate: k_ij = k0 · min(1, exp(-β·ΔH)). Per-species k0 (sets physical time unit; user sets it from D or attempt-frequency).                                  
  - Move encoding for BulkMonomer: each move is (speciesIdx, particleId, dirIdx ∈ 1..6). Slot ID = stable function of these. 6 slots per particle,             
  pre-allocated.                                                                                                                                               
  - Validity: a move with destination occupied (any species) gets rate 0 — kept in the slot, not removed. Move with rate 0 is just absent from the cumulative  
  sum.                                                                                                                                                         
  - Data structures:
    - rates::Vector{Float64} — slot-indexed, length 6·Σ(species.length).                                                                                       
    - fenwick::Vector{Float64} — Fenwick tree over rates (own minimal helpers in gillespie.jl, no new dep).                                                    
    - The lattice itself is the inverted index: given a site x, the species' latS[x] gives the particle at x in O(1). No new index needed.                     
  - Affected-moves rule (single-species BulkMonomer hop p_old → p_new):                                                                                        
    - dirty_sites = NN1(p_old) ∪ NN1(p_new) (≤14 sites).                                                                                                       
    - For each occupied x ∈ dirty_sites: refresh that particle's 6 moves.                                                                                      
    - For each empty x ∈ dirty_sites: for each y ∈ NN(x), if y is occupied and y ∉ dirty_sites, refresh the one move (particle@y, direction y→x).              
    - Multi-species: same, but iterate over all species at each site.                                                                                          
  - Numerical hygiene: full Q rebuild every N_rebuild events (configurable; default e.g. 10⁶) to bound drift.                                                  
  - Sampling: at fixed Δt_sample of simulation time. Wall-clock cap from target! unchanged.                                                                    
                                                                                                                                                               
  ---                                                                                                                                                          
  File layout                                                                                                                                                  
                  
  - src/gillespie.jl (already a stub) — owns KMCSystem, build/update/sample/apply, Fenwick helpers, KMC loop primitives.
  - src/types.jl — add KMCSystem struct and KMCMove if needed (kept minimal; prefer SoA over a struct per move).                                               
  - test/test_kmc_*.jl — one file per phase.                                                                                                                   
  - spec/kmc_diffusion.jl — minimal v1 spec (a single BulkMonomer doing free diffusion; MSD t¹ check).                                                         
  - sim/sweep.jl — untouched; KMC is a different target! body inside the spec, not a new sweep.                                                                
                                                                                                                                                               
  ---                                                                                                                                                          
  Phased plan                                                                                                                                                  
                  
  Phase 1 — scaffolding & rate primitives (no sampling, no time)
                                                                                                                                                               
  Build:                                                                                                                                                       
  - KMCSystem struct: holds slist, rates, fenwick, t_sim::Float64, Q::Float64, beta::Float64, plus per-species k0 (Vector keyed by species index).             
  - slotId(sysIdx, particleId, dir) — pure, inlinable.                                                                                                         
  - computeRate!(sys, slot) — for the (species, particle, dir) decoded from slot: check destination free across all species in slist; if blocked, rate=0; else
  compute ΔH using existing deltaHam!-style logic (we'll likely need a thin local deltaH_kmc helper since deltaHam! is tightly coupled to Metropolis state — to
   be decided: reuse vs minimal duplication; aim for reuse).                                                                                                   
  - buildMoves!(sys) — full sweep, fills rates, builds Fenwick.                                                                                                
                                                                                                                                                               
  Tests (test_kmc_phase1.jl):                                                                                                                                  
  - A 4³ lattice with 1 BulkMonomer, 1 particle, no neighbors → all 6 rates = k0·1.                                                                            
  - 1 particle with 1 occupied NN, intPot=ε, β=1 → moves toward occupant blocked (rate 0), moves away have rate k0·min(1, exp(+βε)) (gain neighbor lost: ΔH=+ε 
  if intPot is binding=negative; sanity-check sign with closed form).                                                                                          
  - Total Q matches sum of rates.                                                                                                                              
                  
  Gate: pass before Phase 2.                                                                                                                                   
                  
  ---                                                                                                                                                          
  Phase 2 — Fenwick sampling
                                                                                                                                                               
  Build:
  - Minimal Fenwick-tree helpers in gillespie.jl: fenwick_update!, fenwick_total, fenwick_find(u) (returns slot whose cumulative sum first exceeds u).         
  - sampleMove(sys) — draws u ~ U(0, Q), returns slot.                                                                                                         
                                                      
  Tests (test_kmc_phase2.jl):                                                                                                                                  
  - Construct rates [1,2,3,4], sample 10⁵ times with fixed seed, χ² test against expected 0.1/0.2/0.3/0.4.                                                     
  - Edge cases: single nonzero rate → always picks it; rate=0 slots never picked.                                                                              
  - fenwick_total ≈ sum(rates) after random updates.                                                                                                           
                                                                                                                                                               
  Gate: pass before Phase 3.                                                                                                                                   
                                                                                                                                                               
  ---                                                                                                                                                          
  Phase 3 — incremental updates
                                                                                                                                                               
  Build:
  - applyMove!(sys, slot) — translate particle, update lattice, update t_sim += -ln(rand())/Q before the rate refresh (Gillespie order).                       
  - refreshAffected!(sys, p_old, p_new) — implements the dirty-site rule above. Each slot whose rate changes: subtract old, recompute, add new,                
  fenwick_update!.
  - Periodic rebuildQ!(sys) for drift correction.                                                                                                              
                  
  Tests (test_kmc_phase3.jl):                                                                                                                                  
  - After every event for 1000 events, compare incremental Q to a from-scratch buildMoves! recomputation. They must agree to ≤1e-12 relative.
  - A "no-op" move (proposed but ΔH refresh is trivial because no neighbors): no rates besides the mover's own should change.                                  
  - Two-particle test: move particle A; verify particle B's 6 rates are correct vs a hand-computed reference for select configs.
                                                                                                                                                               
  Gate: pass before Phase 4. This is the highest-risk phase — most KMC bugs live here.                                                                         
                                                                                                                                                               
  ---                                                                                                                                                          
  Phase 4 — main loop + sim-time sampling                                                                                                                      
                                                                                                                                                               
  Build:
  - kmc_step!(sys) — sample → advance time → applyMove! → refresh.                                                                                             
  - kmc_run!(sys; t_max_wall, t_burn_in, dt_sample, sample_callback):                                                                                          
    - Outer loop bounded by wall-clock (time_ns()).                  
    - Internal t_next_sample accumulator; whenever t_sim ≥ t_next_sample, invoke callback and bump t_next_sample += dt_sample.                                 
    - Burn-in measured in t_sim, not events.                                                                                                                   
  - Output dict additions: t_sim_final, n_events, dt_sample, k0 per species. Existing fields (runtime_min, boxlength, etc.) preserved.                         
                                                                                                                                                               
  Tests (test_kmc_phase4.jl):                                                                                                                                  
  - Equilibrium parity with Metropolis (the gating test). Single BulkMonomer species, 8³ box, modest density and intPot. Run Metropolis to equilibrium; collect
   histogram of ham. Run KMC with same parameters from the same state; collect same histogram. KS-test or moments agree within MC noise (concrete tolerance:   
  |μ_KMC − μ_M| < 3·max(σ_KMC, σ_M)/√Nbatches).                                                                                                                
  - Free-particle MSD: 1 particle, no interactions, large box. MSD vs t_sim should be linear, slope = 2·D where D = k0 · (1²)/6 for 3D simple-cubic hopping    
  with rate k0/6 per direction (i.e. D = k0/6). Tolerance: 5%.                                                                                             
  - Sim-time advancement sanity: Δt = -ln(u)/Q is correct distribution → mean t_sim per event ≈ 1/⟨Q⟩.                                                         
   
  Gate: parity test must pass.                                                                                                                                 
                  
  ---                                                                                                                                                          
  Phase 5 — multi-species cross-interactions
                                                                                                                                                               
  Once Phase 4 holds, multi-species mostly "just works" because deltaH_kmc already iterates intSpecies. The only new bit is that dirty_sites updates may flip
  rates across species (a TF moving changes rates of an RNAP near it).                                                                                         
                  
  Build:                                                                                                                                                       
  - Generalize buildMoves! and refreshAffected! to iterate over all species in slist for both occupancy checks and rate refresh.
  - Confirm intMat cross-terms enter deltaH_kmc correctly.                                                                                                     
   
  Tests (test_kmc_phase5.jl):                                                                                                                                  
  - Two species, attractive cross-interaction, equilibrium contact-frequency matches Metropolis (same parity protocol as Phase 4).
  - Particle of species A moving toward a B: A's rate to that destination is blocked; A's rate parallel to B is enhanced (binding).                            
                  
  Gate: parity test must pass before writing the spec.                                                                                                         
                  
  ---                                                                                                                                                          
  Phase 6 — spec + driver integration
                                                                                                                                                               
  Build:          
  - spec/kmc_diffusion.jl: copied from the closest existing BulkMonomer spec (likely bulk_diffusion_search.jl), target! body swapped to call kmc_run!. Mirror
  the existing output dict shape (species_name, nframes, boxlength, runtime_min, plus the new KMC fields).                                                     
  - Initial spec runs t_max = 1 minute, ncpus = 1 per CRITICAL_PRINCIPLES.                                
                                                                                                                                                               
  Tests:                                                                                                                                                       
  - The spec executes end-to-end via runJL sim/sweep.jl spec/kmc_diffusion.jl 1 data/.../kmc_diffusion_test.
  - Output JSON loads, has expected keys, MSD-from-coords matches the in-memory MSD from Phase 4.                                                              
                                                                                                 
  Gate: actual runJL run completes and produces non-empty data.                                                                                                
                                                                                                                                                               
  ---                                                                                                                                                          
  Phase 7 — global consistency                                                                                                                                 
                              
  - Run the full test suite once (per style_guide.md, only at the end).
  - Add test_kmc_*.jl includes to runtests.jl.                                                                                                                 
  - Confirm no regression in existing Metropolis tests.                                                                                                        
  - If Polymer KMC (v2) is the next step, this is the natural commit/branch boundary.                                                                          
                                                                                                                                                               
  - Free-particle MSD: 1 particle, no interactions, large box. MSD vs t_sim should be linear, slope = 2·D where D = k0 · (1²)/6
  for 3D simple-cubic hopping with rate k0/6 per direction (i.e. D = k0/6). Tolerance: 5%.
  - Sim-time advancement sanity: Δt = -ln(u)/Q is correct distribution → mean t_sim per event ≈ 1/⟨Q⟩.

  Gate: parity test must pass.

  ---
  Phase 5 — multi-species cross-interactions

  Once Phase 4 holds, multi-species mostly "just works" because deltaH_kmc already iterates intSpecies. The only new bit is that
  dirty_sites updates may flip rates across species (a TF moving changes rates of an RNAP near it).

  Build:
  - Generalize buildMoves! and refreshAffected! to iterate over all species in slist for both occupancy checks and rate refresh.
  - Confirm intMat cross-terms enter deltaH_kmc correctly.

  Tests (test_kmc_phase5.jl):
  - Two species, attractive cross-interaction, equilibrium contact-frequency matches Metropolis (same parity protocol as Phase 4).
  - Particle of species A moving toward a B: A's rate to that destination is blocked; A's rate parallel to B is enhanced
  (binding).

  Gate: parity test must pass before writing the spec.

  ---
  Phase 6 — spec + driver integration

  Build:
  - spec/kmc_diffusion.jl: copied from the closest existing BulkMonomer spec (likely bulk_diffusion_search.jl), target! body
  Build:
  - spec/kmc_diffusion.jl: copied from the closest existing BulkMonomer spec (likely bulk_diffusion_search.jl), target! body swapped to call kmc_run!. Mirror
  the existing output dict shape (species_name, nframes, boxlength, runtime_min, plus the new KMC fields).
  - Initial spec runs t_max = 1 minute, ncpus = 1 per CRITICAL_PRINCIPLES.

  Tests:
  - The spec executes end-to-end via runJL sim/sweep.jl spec/kmc_diffusion.jl 1 data/.../kmc_diffusion_test.
  - Output JSON loads, has expected keys, MSD-from-coords matches the in-memory MSD from Phase 4.

  Gate: actual runJL run completes and produces non-empty data.

  ---
  Phase 7 — global consistency

  - Run the full test suite once (per style_guide.md, only at the end).
  - Add test_kmc_*.jl includes to runtests.jl.
  - Confirm no regression in existing Metropolis tests.
  - If Polymer KMC (v2) is the next step, this is the natural commit/branch boundary.