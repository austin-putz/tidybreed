# Database-side helpers for the term/member genome-effect model.
#
# helper-genome-effects.R is Phase A: pure R, no database, the hand-computed
# origin truth table. This file is the other side — small utilities for tests
# that go through the real tables.

#' A flat view of the generated additive effects
#'
#' Several suites want the old one-row-per-(trait, locus, line) shape to assert
#' against: it is the natural grain for "did this locus get this coefficient and
#' this centre?". Reconstructing it from the three tables in a view keeps those
#' assertions readable **and** keeps them honest — the view is derived from the
#' real storage, so a change in the storage breaks it rather than silently
#' passing. Restricted to the reserved owner, which is the only owner whose
#' terms are guaranteed order-one additive.
ge_flat_view <- function(pop) {
  DBI::dbExecute(pop$db_conn, paste0(
    "CREATE OR REPLACE VIEW gen_add_flat AS ",
    "SELECT e.trait_name, l.locus_name, m.locus_id, ",
    "       m.center_value AS base_allele_freq, e.genome_value, ",
    "       o.line_name, o.parent_origin ",
    "FROM genome_effects e ",
    "JOIN genome_effect_members m USING (id_genome_effect) ",
    "JOIN genome_effect_loci   l USING (id_genome_effect, member_slot) ",
    "LEFT JOIN genome_effect_member_origins o ",
    "  USING (id_genome_effect, member_slot) ",
    "WHERE e.effect_owner = 'generated_additive_tbv'"))
  pop
}
