//! Decision-oriented primitives for lipid-aware lead-optimisation campaigns.
//!
//! This module deliberately separates campaign analysis from the simulation
//! engine. FerrumMD (or an external engine such as OpenMM) can generate the
//! measurements, while these types rank membrane compositions, R-group
//! proposals, and fragment molecular orbital (FMO) interactions consistently.

use std::cmp::Ordering;
use std::collections::BTreeMap;

/// A bilayer composition used for a simulation or experimental condition.
#[derive(Clone, Debug, PartialEq)]
pub struct MembraneEnvironment {
    pub name: String,
    /// Lipid name to mole fraction. Fractions are expected to sum to one.
    pub lipid_fractions: BTreeMap<String, f64>,
}

impl MembraneEnvironment {
    pub fn new(
        name: impl Into<String>,
        lipid_fractions: impl IntoIterator<Item = (impl Into<String>, f64)>,
    ) -> Result<Self, LeadOptimizationError> {
        let lipid_fractions = lipid_fractions
            .into_iter()
            .map(|(lipid, fraction)| (lipid.into(), fraction))
            .collect::<BTreeMap<_, _>>();

        if lipid_fractions.is_empty()
            || lipid_fractions
                .values()
                .any(|fraction| !fraction.is_finite() || *fraction < 0.0)
        {
            return Err(LeadOptimizationError::InvalidMembraneComposition);
        }

        let total: f64 = lipid_fractions.values().sum();
        if (total - 1.0).abs() > 1e-6 {
            return Err(LeadOptimizationError::InvalidMembraneComposition);
        }

        Ok(Self {
            name: name.into(),
            lipid_fractions,
        })
    }
}

#[derive(Clone, Debug, PartialEq)]
pub enum LeadOptimizationError {
    InvalidMembraneComposition,
    InvalidWeight,
    EmptyCampaign,
}

/// Simulation measurements for one compound in one membrane environment.
#[derive(Clone, Debug, PartialEq)]
pub struct MembranePermutationResult {
    pub compound_id: String,
    pub membrane: MembraneEnvironment,
    /// Membrane-crossing free-energy barrier in kJ/mol; lower is preferred.
    pub permeation_barrier: f64,
    /// Binding free energy in kJ/mol; more negative is preferred.
    pub binding_free_energy: f64,
    /// Bilayer disruption or instability penalty; lower is preferred.
    pub membrane_disruption: f64,
    /// Statistical uncertainty associated with the result.
    pub uncertainty: f64,
}

/// Weights for translating simulation outputs into a campaign ranking.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ObjectiveWeights {
    pub permeation: f64,
    pub binding: f64,
    pub disruption: f64,
    pub uncertainty: f64,
}

impl Default for ObjectiveWeights {
    fn default() -> Self {
        Self {
            permeation: 1.0,
            binding: 1.0,
            disruption: 1.0,
            uncertainty: 0.5,
        }
    }
}

impl ObjectiveWeights {
    fn validate(self) -> Result<Self, LeadOptimizationError> {
        let values = [
            self.permeation,
            self.binding,
            self.disruption,
            self.uncertainty,
        ];
        if values
            .iter()
            .any(|value| !value.is_finite() || *value < 0.0)
        {
            Err(LeadOptimizationError::InvalidWeight)
        } else {
            Ok(self)
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct RankedMembraneResult {
    pub result: MembranePermutationResult,
    pub score: f64,
}

/// Rank membrane permutations without hiding uncertainty or bilayer damage.
pub fn rank_membrane_permutations(
    results: Vec<MembranePermutationResult>,
    weights: ObjectiveWeights,
) -> Result<Vec<RankedMembraneResult>, LeadOptimizationError> {
    if results.is_empty() {
        return Err(LeadOptimizationError::EmptyCampaign);
    }
    let weights = weights.validate()?;

    let mut ranked = results
        .into_iter()
        .map(|result| RankedMembraneResult {
            score: weights.permeation * result.permeation_barrier
                + weights.binding * result.binding_free_energy
                + weights.disruption * result.membrane_disruption
                + weights.uncertainty * result.uncertainty,
            result,
        })
        .collect::<Vec<_>>();
    ranked.sort_by(|left, right| finite_cmp(left.score, right.score));
    Ok(ranked)
}

/// FEgrow-compatible proposal augmented with membrane-aware observables.
#[derive(Clone, Debug, PartialEq)]
pub struct RGroupProposal {
    pub compound_id: String,
    pub attachment_atom: usize,
    pub r_group_smiles: String,
    pub pocket_interaction_energy: f64,
    pub membrane_transfer_penalty: f64,
    pub steric_clash_penalty: f64,
    pub synthetic_accessibility_penalty: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct RankedRGroupProposal {
    pub proposal: RGroupProposal,
    pub score: f64,
}

/// Rank R-group growth proposals using both pocket and bilayer terms.
///
/// Lower scores are better. This provides the lipid-aware scoring seam for a
/// FEgrow adapter while leaving conformer generation to FEgrow itself.
pub fn rank_r_group_proposals(
    proposals: Vec<RGroupProposal>,
) -> Result<Vec<RankedRGroupProposal>, LeadOptimizationError> {
    if proposals.is_empty() {
        return Err(LeadOptimizationError::EmptyCampaign);
    }

    let mut ranked = proposals
        .into_iter()
        .map(|proposal| RankedRGroupProposal {
            score: proposal.pocket_interaction_energy
                + proposal.membrane_transfer_penalty
                + proposal.steric_clash_penalty
                + proposal.synthetic_accessibility_penalty,
            proposal,
        })
        .collect::<Vec<_>>();
    ranked.sort_by(|left, right| finite_cmp(left.score, right.score));
    Ok(ranked)
}

#[derive(Clone, Copy, Debug, Eq, PartialEq, Ord, PartialOrd)]
pub enum FragmentKind {
    Ligand,
    Protein,
    Lipid,
    Water,
    Cofactor,
}

/// Pair interaction energy from an external FMO calculation.
#[derive(Clone, Debug, PartialEq)]
pub struct FmoPairInteraction {
    pub fragment_a: String,
    pub kind_a: FragmentKind,
    pub fragment_b: String,
    pub kind_b: FragmentKind,
    /// Pair interaction energy in kJ/mol; negative values are favourable.
    pub interaction_energy: f64,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct FmoInteractionSummary {
    pub protein: f64,
    pub lipid: f64,
    pub water: f64,
    pub cofactor: f64,
    pub total_environment: f64,
    pub strongest_favourable_contacts: Vec<FmoPairInteraction>,
}

/// Summarise ligand-environment FMO terms, explicitly retaining lipid contacts.
pub fn summarize_ligand_fmo(
    interactions: &[FmoPairInteraction],
    strongest_contact_count: usize,
) -> FmoInteractionSummary {
    let mut summary = FmoInteractionSummary::default();
    let mut ligand_contacts = Vec::new();

    for interaction in interactions {
        let environment_kind = match (interaction.kind_a, interaction.kind_b) {
            (FragmentKind::Ligand, other) if other != FragmentKind::Ligand => Some(other),
            (other, FragmentKind::Ligand) if other != FragmentKind::Ligand => Some(other),
            _ => None,
        };

        if let Some(kind) = environment_kind {
            match kind {
                FragmentKind::Protein => summary.protein += interaction.interaction_energy,
                FragmentKind::Lipid => summary.lipid += interaction.interaction_energy,
                FragmentKind::Water => summary.water += interaction.interaction_energy,
                FragmentKind::Cofactor => summary.cofactor += interaction.interaction_energy,
                FragmentKind::Ligand => {}
            }
            summary.total_environment += interaction.interaction_energy;
            ligand_contacts.push(interaction.clone());
        }
    }

    ligand_contacts
        .sort_by(|left, right| finite_cmp(left.interaction_energy, right.interaction_energy));
    ligand_contacts.truncate(strongest_contact_count);
    summary.strongest_favourable_contacts = ligand_contacts;
    summary
}

fn finite_cmp(left: f64, right: f64) -> Ordering {
    left.partial_cmp(&right).unwrap_or(Ordering::Equal)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn membrane(name: &str) -> MembraneEnvironment {
        MembraneEnvironment::new(name, [("POPC", 0.7), ("CHOL", 0.3)]).unwrap()
    }

    #[test]
    fn rejects_non_normalized_membranes() {
        assert_eq!(
            MembraneEnvironment::new("invalid", [("POPC", 0.8)]),
            Err(LeadOptimizationError::InvalidMembraneComposition)
        );
    }

    #[test]
    fn membrane_ranking_balances_binding_and_permeation() {
        let ranked = rank_membrane_permutations(
            vec![
                MembranePermutationResult {
                    compound_id: "CMP-1".into(),
                    membrane: membrane("reference"),
                    permeation_barrier: 25.0,
                    binding_free_energy: -40.0,
                    membrane_disruption: 2.0,
                    uncertainty: 1.0,
                },
                MembranePermutationResult {
                    compound_id: "CMP-2".into(),
                    membrane: membrane("reference"),
                    permeation_barrier: 18.0,
                    binding_free_energy: -38.0,
                    membrane_disruption: 1.0,
                    uncertainty: 0.5,
                },
            ],
            ObjectiveWeights::default(),
        )
        .unwrap();

        assert_eq!(ranked[0].result.compound_id, "CMP-2");
    }

    #[test]
    fn r_group_ranking_penalizes_a_lipid_incompatible_design() {
        let ranked = rank_r_group_proposals(vec![
            RGroupProposal {
                compound_id: "lipophilic".into(),
                attachment_atom: 4,
                r_group_smiles: "c1ccccc1".into(),
                pocket_interaction_energy: -15.0,
                membrane_transfer_penalty: 12.0,
                steric_clash_penalty: 0.0,
                synthetic_accessibility_penalty: 1.0,
            },
            RGroupProposal {
                compound_id: "balanced".into(),
                attachment_atom: 4,
                r_group_smiles: "CO".into(),
                pocket_interaction_energy: -10.0,
                membrane_transfer_penalty: 2.0,
                steric_clash_penalty: 0.0,
                synthetic_accessibility_penalty: 1.0,
            },
        ])
        .unwrap();

        assert_eq!(ranked[0].proposal.compound_id, "balanced");
    }

    #[test]
    fn fmo_summary_separates_lipid_from_protein_contacts() {
        let summary = summarize_ligand_fmo(
            &[
                FmoPairInteraction {
                    fragment_a: "LIG".into(),
                    kind_a: FragmentKind::Ligand,
                    fragment_b: "ASP42".into(),
                    kind_b: FragmentKind::Protein,
                    interaction_energy: -20.0,
                },
                FmoPairInteraction {
                    fragment_a: "CHOL7".into(),
                    kind_a: FragmentKind::Lipid,
                    fragment_b: "LIG".into(),
                    kind_b: FragmentKind::Ligand,
                    interaction_energy: -6.0,
                },
            ],
            1,
        );

        assert_eq!(summary.protein, -20.0);
        assert_eq!(summary.lipid, -6.0);
        assert_eq!(summary.total_environment, -26.0);
        assert_eq!(summary.strongest_favourable_contacts[0].fragment_b, "ASP42");
    }
}
