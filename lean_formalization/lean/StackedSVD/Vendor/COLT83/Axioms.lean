/-
Copyright (c) 2026 Rémy Degenne. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Rémy Degenne
-/
-- Backported from RemyDegenne/colt-2026-83, commit 2d846025a86ee6760a805901c04237039f987457
-- (Mathlib v4.34.0-rc2), to Mathlib v4.33.0 on 2026-08-29.
module

public import StackedSVD.Vendor.COLT83.Mathlib.Probability.SteinReal
public import StackedSVD.Vendor.COLT83.Mathlib.Probability.SteinIdentity
public import StackedSVD.Vendor.COLT83.Mathlib.Probability.GaussianInterpolation
public import StackedSVD.Vendor.COLT83.Mathlib.Probability.SudakovFernique
public import StackedSVD.Vendor.COLT83.Mathlib.Probability.BorellTIS

set_option autoImplicit false

/-!
# Axiom audit of the backported COLT83 Gaussian files

This file imports the five headline modules and prints the axioms of the two results the
project depends on. The expected output for each is `propext, Classical.choice, Quot.sound`.
-/

#print axioms ProbabilityTheory.sudakov_fernique
#print axioms ProbabilityTheory.integral_inner_mul_stdGaussian
