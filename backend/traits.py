"""
traits.py -- ABO/RhD blood group prediction and well-established trait calls.

Two independent pieces live here:

  1. Blood group prediction from array-accessible tag SNPs. Consumer arrays do
     not genotype ABO cleanly, so every function here reports what it could not
     resolve as loudly as what it could. A confident call is never returned
     when the decisive tag is missing.
  2. A declarative TRAITS table covering only non-medical, well-replicated,
     array-accessible traits, plus a small total evaluator for it. Traits are
     neutral: nothing in this module is coloured good or bad.

Nothing here is diagnostic. Blood group in particular must never be used for
transfusion or pregnancy decisions; serology is the only valid source for that.
"""

from typing import Any

try:  # orientation.py may be added to the package after this module
    from backend import orientation as _orientation
except ImportError:  # pragma: no cover - only before orientation.py lands
    try:
        import orientation as _orientation  # type: ignore
    except ImportError:
        _orientation = None  # type: ignore


_NOCALL_ALLELES = {"", "N", "-", "--", "0", "00", "?", "."}
NOCALL_KEY = "NN"

_CONFIDENCE_ORDER = ("none", "low", "moderate", "high")


def _sorted_pair(a1: Any, a2: Any) -> tuple[str, str]:
    """Return the two alleles in canonical order, via orientation when available."""
    left = str(a1 if a1 is not None else "").strip().upper()
    right = str(a2 if a2 is not None else "").strip().upper()
    if _orientation is not None:
        try:
            pair = _orientation.sort_alleles(left, right)
            if isinstance(pair, (list, tuple)) and len(pair) == 2:
                first, second = str(pair[0]).upper(), str(pair[1]).upper()
                if sorted((first, second)) == sorted((left, right)):
                    return first, second
        except Exception:
            pass
    ordered = sorted((left, right))
    return ordered[0], ordered[1]


def _key_from_text(text: Any) -> str:
    """Normalise a 2-character genotype string to a sorted key."""
    raw = str(text if text is not None else "").strip().upper()
    if len(raw) != 2:
        return NOCALL_KEY
    first, second = _sorted_pair(raw[0], raw[1])
    return f"{first}{second}"


def _index(genotypes: Any) -> dict:
    """Lower-case the rsID keys of the caller's genotype mapping."""
    if not isinstance(genotypes, dict):
        return {}
    return {str(k).strip().lower(): v for k, v in genotypes.items()}


def _read(genotypes: dict, rsid: str) -> str | None:
    """Return the canonical genotype key for one rsID, or None when absent.

    A genotyped but uncalled position returns "NN"; a position the array does
    not carry at all returns None. The distinction matters because a missing
    tag and a failed tag need different wording in the caveat.
    """
    if rsid not in genotypes:
        return None
    value = genotypes[rsid]
    if isinstance(value, (list, tuple)) and len(value) == 2:
        first, second = str(value[0]).upper().strip(), str(value[1]).upper().strip()
    else:
        raw = str(value if value is not None else "").strip().upper()
        if len(raw) != 2:
            return NOCALL_KEY
        first, second = raw[0], raw[1]
    if first in _NOCALL_ALLELES or second in _NOCALL_ALLELES:
        return NOCALL_KEY
    ordered = _sorted_pair(first, second)
    return f"{ordered[0]}{ordered[1]}"


def _has_call(genotypes: dict, rsid: str) -> bool:
    """True when the rsID is present and is a real call."""
    key = _read(genotypes, rsid)
    return key is not None and key != NOCALL_KEY


def _copies(genotypes: dict, rsid: str, allele: str) -> int | None:
    """Count copies (0, 1 or 2) of one allele, or None when there is no call."""
    key = _read(genotypes, rsid)
    if key is None or key == NOCALL_KEY:
        return None
    return sum(1 for base in key if base == str(allele).upper())


def _weaker(first: str, second: str) -> str:
    """Return the weaker of two confidence labels."""
    try:
        return _CONFIDENCE_ORDER[min(
            _CONFIDENCE_ORDER.index(first), _CONFIDENCE_ORDER.index(second)
        )]
    except ValueError:
        return "none"


def _downgrade(confidence: str, steps: int = 1) -> str:
    """Move a confidence label down the ladder, never below "none"."""
    try:
        index = _CONFIDENCE_ORDER.index(confidence)
    except ValueError:
        return "none"
    return _CONFIDENCE_ORDER[max(0, index - steps)]


# ---------------------------------------------------------------------------
# ABO blood group
# ---------------------------------------------------------------------------
# rs8176719  exon 6 single-base deletion. A deletion on both copies abolishes
#            transferase activity and gives group O. This is the ONLY decisive
#            tag; without it no group can be called.
# rs8176746  Leu266Met, distinguishes the A from the B transferase.
# rs8176747  Gly268Ala, second A versus B discriminator.
# rs505922   intron 1, in strong linkage disequilibrium with the O allele.
#            Used only as a consistency check, never to make the call.
# rs1053878  A1 versus A2 subgroup. Rarely on consumer arrays, so subgroups
#            are reported as unresolved rather than guessed.
ABO_DECISIVE_TAG = "rs8176719"
ABO_AB_TAGS = ("rs8176746", "rs8176747")
ABO_LD_TAG = "rs505922"
ABO_SUBGROUP_TAG = "rs1053878"
ABO_TAGS = (ABO_DECISIVE_TAG, "rs8176746", "rs8176747", ABO_LD_TAG)

# Deletion (O) and intact (A/B backbone) tokens seen at rs8176719 across
# vendors. Several arrays emit D/I rather than bases for this indel.
#
# "-" is deliberately NOT a deletion token. 23andMe and AncestryDNA both write
# "--" for a probe that failed, and a single "-" for a half call. Treating "-"
# as the exon 6 deletion turns the most common failure mode of the one decisive
# ABO tag into a confidently reported group O, which is exactly the false
# positive this module is supposed to refuse to produce. Vendors that genuinely
# report this indel use D and I.
_DELETION_TOKENS = {"D", "DEL"}
_INTACT_TOKENS = {"I", "G", "C", "INS"}

# Tokens that mean "no result here" at an indel position. Checked before the
# deletion and intact sets so a failed probe can never be read as a genotype.
_INDEL_NOCALL_TOKENS = {"", "-", "--", "0", "00", "N", "NN", "?", "??"}

# Allele that rides with the B transferase at each A/B discriminator, on the
# strand the major vendors report.
_B_ALLELES = {"rs8176746": "T", "rs8176747": "C"}


def _read_indel(genotypes: dict, rsid: str) -> tuple[int, str] | None:
    """Read rs8176719 style indel calls as (deletion_copies, raw_token).

    Handles the D/I token style and the intact-base style (G or C). A "-" is
    treated as a NO-CALL, not as the deletion allele, because that is what the
    major vendors mean by it. Returns None when the position is absent,
    uncalled, or written with tokens this function cannot interpret, which is a
    very common outcome for this indel and is reported honestly rather than
    guessed.
    """
    if rsid not in genotypes:
        return None
    value = genotypes[rsid]
    if isinstance(value, (list, tuple)) and len(value) == 2:
        alleles = [str(value[0]).strip().upper(), str(value[1]).strip().upper()]
    else:
        raw = str(value if value is not None else "").strip().upper()
        # Reject whole-token no-calls first. "--" would otherwise split into
        # two "-" alleles and, before this guard existed, read as a homozygous
        # deletion and produce a confident group O from a failed probe.
        if raw in _INDEL_NOCALL_TOKENS:
            return None
        if len(raw) != 2:
            return None
        alleles = [raw[0], raw[1]]
    # Per-allele no-call guard, so a half call such as ("D", "-") is also
    # refused rather than counted as one deletion copy.
    if any(a in _INDEL_NOCALL_TOKENS for a in alleles):
        return None
    deletions = 0
    for allele in alleles:
        if allele in _DELETION_TOKENS:
            deletions += 1
        elif allele not in _INTACT_TOKENS:
            return None
    return deletions, "".join(sorted(alleles))


def _abo_b_copies(genotypes: dict) -> tuple[int | None, list[str], bool]:
    """Count B-transferase allele copies from the A/B discriminators.

    Returns (copies, tags_used, consistent). ``copies`` is None when neither
    discriminator was called. ``consistent`` is False when the two tags
    disagree, which normally means a strand mismatch rather than a real
    biological oddity.
    """
    counts: dict[str, int] = {}
    for rsid in ABO_AB_TAGS:
        copies = _copies(genotypes, rsid, _B_ALLELES[rsid])
        if copies is not None:
            counts[rsid] = copies
    if not counts:
        return None, [], True
    values = set(counts.values())
    consistent = len(values) == 1
    return max(counts.values()), sorted(counts.keys()), consistent


_ABO_GROUP_OF = {"OO": "O", "AO": "A", "AA": "A", "BO": "B", "BB": "B", "AB": "AB"}


def predict_abo(genotypes: dict) -> dict:
    """Predict ABO blood group from array-accessible tag SNPs.

    ``genotypes`` maps lower-case rsID to either a 2-item sequence of alleles
    or a 2-character genotype string.

    Honest limitations, all of which are reported back to the caller:
    consumer arrays do not genotype ABO cleanly; rs8176719 is a single-base
    deletion that many arrays report as a no-call or as D/I tokens rather than
    as bases; and the A1 versus A2 subgroup (rs1053878) is usually not
    resolvable at all. When rs8176719 is missing or uncalled the result is
    always "unknown" with confidence "none", never a guess.
    """
    g = _index(genotypes)
    indel = _read_indel(g, ABO_DECISIVE_TAG)
    missing: list[str] = []
    if indel is None:
        missing.append(ABO_DECISIVE_TAG)
    for rsid in ABO_AB_TAGS + (ABO_LD_TAG, ABO_SUBGROUP_TAG):
        if not _has_call(g, rsid):
            missing.append(rsid)

    b_copies, tags_used, consistent = _abo_b_copies(g)
    deciding: list[dict] = []
    alternatives: list[dict] = []

    if indel is None:
        reason = (
            f"{ABO_DECISIVE_TAG} is the only tag that identifies the O allele, and it "
            "is a single-base deletion rather than a substitution. Consumer arrays "
            "commonly omit it, return a no-call, or write it as D/I tokens that cannot "
            "be read as bases, so no blood group can be assigned."
        )
        if b_copies is None:
            alternatives = [
                {"abo": group, "why": "No ABO tag resolved, so every group remains possible."}
                for group in ("O", "A", "B", "AB")
            ]
        elif b_copies == 0:
            alternatives = [
                {"abo": "A", "why": "No B allele seen at the A/B tags, so A or O remain."},
                {"abo": "O", "why": "O cannot be excluded without rs8176719."},
            ]
        elif b_copies == 1:
            alternatives = [
                {"abo": "AB", "why": "One B allele seen; the other copy may be A."},
                {"abo": "B", "why": "One B allele seen; the other copy may be O."},
            ]
        else:
            alternatives = [
                {"abo": "B", "why": "Two B alleles at the A/B tags."},
                {"abo": "O", "why": "O cannot be excluded without rs8176719."},
            ]
        return {
            "abo": "unknown",
            "genotype_call": "",
            "confidence": "none",
            "deciding": deciding,
            "missing": missing,
            "caveat": reason,
            "alternatives": alternatives,
        }

    del_copies, del_token = indel
    deciding.append({
        "rsid": ABO_DECISIVE_TAG,
        "genotype": del_token,
        "role": f"O allele tag, exon 6 deletion, {del_copies} deleted copy/copies",
    })
    for rsid in tags_used:
        deciding.append({
            "rsid": rsid,
            "genotype": _read(g, rsid) or "",
            "role": f"A versus B discriminator, B allele {_B_ALLELES[rsid]}",
        })

    caveats: list[str] = []
    if del_copies == 2:
        genotype_call, abo, confidence = "OO", "O", "high"
        caveats.append(
            "The exon 6 deletion is present on both copies, which abolishes ABO "
            "transferase activity, so group O follows directly and the A/B tags are "
            "not needed."
        )
    elif b_copies is None:
        genotype_call, abo, confidence = "", "unknown", "low"
        caveats.append(
            f"{ABO_DECISIVE_TAG} was read, so at least one intact ABO allele is present, "
            f"but neither {ABO_AB_TAGS[0]} nor {ABO_AB_TAGS[1]} was called, so A cannot be "
            "distinguished from B."
        )
        if del_copies == 1:
            alternatives = [
                {"abo": "A", "why": "Genotype AO if the intact allele is an A transferase."},
                {"abo": "B", "why": "Genotype BO if the intact allele is a B transferase."},
            ]
        else:
            alternatives = [
                {"abo": "A", "why": "Genotype AA if both intact alleles are A."},
                {"abo": "AB", "why": "Genotype AB if one allele is A and one is B."},
                {"abo": "B", "why": "Genotype BB if both intact alleles are B."},
            ]
    else:
        if del_copies == 1:
            genotype_call = "BO" if b_copies >= 1 else "AO"
        elif b_copies == 0:
            genotype_call = "AA"
        elif b_copies == 1:
            genotype_call = "AB"
        else:
            genotype_call = "BB"
        abo = _ABO_GROUP_OF[genotype_call]
        confidence = "high" if len(tags_used) == 2 and consistent else "moderate"
        if not consistent:
            confidence = "low"
            caveats.append(
                f"{ABO_AB_TAGS[0]} and {ABO_AB_TAGS[1]} disagree on the number of B alleles, "
                "which usually means one of them was reported on the opposite strand."
            )
            alternatives = [
                {"abo": "A", "why": "Consistent if the A/B tags are read on one strand."},
                {"abo": "B", "why": "Consistent if the A/B tags are read on the other strand."},
            ]
        if len(tags_used) == 1:
            caveats.append(
                f"Only {tags_used[0]} was available to separate A from B; the second "
                "discriminator would normally confirm it."
            )

    ld_copies = _copies(g, ABO_LD_TAG, "C")
    if ld_copies is None:
        caveats.append(
            f"{ABO_LD_TAG}, the intron 1 tag in strong linkage with the O allele, was not "
            "available as a consistency check."
        )
    else:
        deciding.append({
            "rsid": ABO_LD_TAG,
            "genotype": _read(g, ABO_LD_TAG) or "",
            "role": "linkage check only, C allele travels with the O allele",
        })
        inconsistent = ld_copies < del_copies or (del_copies == 0 and ld_copies == 2)
        if inconsistent:
            confidence = _downgrade(confidence)
            caveats.append(
                f"{ABO_LD_TAG} carries {ld_copies} C allele(s), which does not fit "
                f"{del_copies} deleted copy/copies at {ABO_DECISIVE_TAG}; linkage between the "
                "two varies by ancestry, so treat the call with extra caution."
            )

    if abo in ("A", "AB"):
        caveats.append(
            f"The A1 versus A2 subgroup depends on {ABO_SUBGROUP_TAG}, which consumer arrays "
            "usually do not carry, so the subgroup is left unresolved."
        )
    caveats.append(
        "Array-based ABO prediction is not a substitute for serological typing and must "
        "never be used for transfusion decisions."
    )

    return {
        "abo": abo,
        "genotype_call": genotype_call,
        "confidence": confidence,
        "deciding": deciding,
        "missing": missing,
        "caveat": " ".join(caveats),
        "alternatives": alternatives,
    }


# ---------------------------------------------------------------------------
# RhD status
# ---------------------------------------------------------------------------
# rs590787 is the commonly used array proxy for the whole-gene RHD deletion in
# European populations. rs586178 sits on the same haplotype and is used only as
# a weak supporting check. Both allele assignments below are strand dependent:
# if a vendor reports the opposite strand the proxy inverts, which is one more
# reason a negative call here is soft.
RH_PRIMARY_TAG = "rs590787"
RH_SUPPORT_TAG = "rs586178"
_RH_DELETION_ALLELES = {"rs590787": "C", "rs586178": "T"}

RH_CAVEAT_BASE = (
    "The common European RHD whole-gene deletion is the main cause of RhD negativity, "
    "and array proxies tag that one deletion haplotype only. Non-European RhD-negative "
    "haplotypes, such as the RHD pseudogene common in African ancestry and the RHD-CE-Ds "
    "hybrid alleles, are usually missed entirely, so a positive call from this proxy is "
    "considerably more trustworthy than a negative one. This is a proxy for a structural "
    "variant, not a direct read of RHD, and it must never be used for transfusion or "
    "pregnancy decisions."
)


def predict_rh(genotypes: dict) -> dict:
    """Predict RhD status from RHD-region tag SNPs.

    rs590787 is the commonly used array proxy for the RHD whole-gene deletion
    in European populations; rs586178 is treated as a weak supporting tag on
    the same haplotype.

    Returns {"rh", "confidence", "deciding", "caveat"}. Confidence is
    deliberately asymmetric: a "positive" call can be reached from a single
    intact RHD copy and is robust, whereas a "negative" call depends on the
    proxy having captured the right deletion haplotype, so it is never rated
    above "moderate".
    """
    g = _index(genotypes)
    deciding: list[dict] = []
    primary = _copies(g, RH_PRIMARY_TAG, _RH_DELETION_ALLELES[RH_PRIMARY_TAG])
    if primary is None:
        return {
            "rh": "unknown",
            "confidence": "none",
            "deciding": deciding,
            "caveat": (
                f"{RH_PRIMARY_TAG}, the array proxy for the RHD deletion, was not genotyped "
                "or returned a no-call, so RhD status cannot be estimated. " + RH_CAVEAT_BASE
            ),
        }

    deciding.append({
        "rsid": RH_PRIMARY_TAG,
        "genotype": _read(g, RH_PRIMARY_TAG) or "",
        "role": f"RHD deletion proxy, {primary} deletion-tag allele(s)",
    })
    support = _copies(g, RH_SUPPORT_TAG, _RH_DELETION_ALLELES[RH_SUPPORT_TAG])
    if support is not None:
        deciding.append({
            "rsid": RH_SUPPORT_TAG,
            "genotype": _read(g, RH_SUPPORT_TAG) or "",
            "role": "supporting haplotype tag",
        })

    notes: list[str] = []
    if primary == 2:
        rh = "negative"
        agrees = support is not None and support == 2
        confidence = "moderate" if agrees else "low"
        notes.append(
            "Both copies carry the deletion-tag allele, which is consistent with the "
            "European RHD deletion haplotype on both chromosomes."
        )
        if support is None:
            notes.append(f"{RH_SUPPORT_TAG} was unavailable to support the negative call.")
        elif not agrees:
            notes.append(f"{RH_SUPPORT_TAG} does not agree with the primary tag.")
    else:
        rh = "positive"
        agrees = support is not None and support <= 1
        confidence = "high" if agrees else "moderate"
        notes.append(
            "At least one chromosome lacks the deletion-tag allele, so at least one intact "
            "RHD gene is expected, and a single intact copy is enough to be RhD positive."
        )
        if support is None:
            notes.append(f"{RH_SUPPORT_TAG} was unavailable as a supporting check.")

    return {
        "rh": rh,
        "confidence": confidence,
        "deciding": deciding,
        "caveat": " ".join(notes + [RH_CAVEAT_BASE]),
    }


def predict_blood_type(genotypes: dict) -> dict:
    """Combine ABO and RhD into a single blood type with the weaker confidence."""
    abo = predict_abo(genotypes)
    rh = predict_rh(genotypes)
    sign = {"positive": "+", "negative": "-"}.get(rh["rh"], "")
    if abo["abo"] != "unknown" and sign:
        blood_type = f"{abo['abo']}{sign}"
    else:
        blood_type = "unknown"
    confidence = _weaker(abo["confidence"], rh["confidence"])
    if blood_type == "unknown":
        summary = (
            "Blood type could not be predicted from this file: "
            f"ABO is {abo['abo']} (confidence {abo['confidence']}) and RhD is {rh['rh']} "
            f"(confidence {rh['confidence']})."
        )
    else:
        summary = (
            f"Array tags are most consistent with blood type {blood_type}, at {confidence} "
            "confidence, and this is a prediction from proxy markers rather than a blood test."
        )
    return {
        "blood_type": blood_type,
        "abo": abo,
        "rh": rh,
        "confidence": confidence,
        "summary": summary,
    }


# ---------------------------------------------------------------------------
# Trait definitions
# ---------------------------------------------------------------------------
# Only well-established, non-medical, array-accessible traits belong here.
# Rule specs are tiny and declarative:
#   {"rsid": "rs4988235", "genotype": ["CC"]}          exact genotype, order free
#   {"rsid": "rs713598", "min_copies": {"allele": "G", "n": 1}}
#   {"rsid": "rs1805007", "max_copies": {"allele": "T", "n": 0}}
#   {"all": [spec, spec]} / {"any": [spec, spec]}      combinators
# Rules are tried in order and the first one that evaluates True wins.

TRAITS: list[dict] = [
    {
        "key": "lactase_persistence",
        "name": "Lactase persistence",
        "category": "Diet",
        "rsids": ["rs4988235"],
        "rules": [
            {"when": {"rsid": "rs4988235", "genotype": ["TT"]},
             "call": "Lactase persistent",
             "detail": "Two copies of the persistence allele, so lactase production normally continues through adult life."},
            {"when": {"rsid": "rs4988235", "genotype": ["CT"]},
             "call": "Likely lactase persistent",
             "detail": "One persistence allele is usually enough to digest lactose, though tolerance is often lower than in TT."},
            {"when": {"rsid": "rs4988235", "genotype": ["CC"]},
             "call": "Lactase non-persistent",
             "detail": "No copy of the European persistence allele, so lactase activity typically falls after weaning and dairy may cause symptoms."},
        ],
        "default_call": "Undetermined",
        "magnitude": 2.0,
        "evidence": "MCM6 -13910C>T is the best replicated lactase persistence variant, with large effects in European and some South Asian groups.",
        "caveat": "This tag explains lactase persistence in Europeans. Several African and Middle Eastern populations achieve persistence through different variants that this SNP does not capture, so a non-persistent call in those ancestries can be wrong.",
    },
    {
        "key": "alcohol_flush",
        "name": "Alcohol flush reaction",
        "category": "Metabolism",
        "rsids": ["rs671"],
        "rules": [
            {"when": {"rsid": "rs671", "genotype": ["AA"]},
             "call": "Strong alcohol flush reaction",
             "detail": "Two ALDH2*2 copies leave almost no acetaldehyde clearance, so even small amounts of alcohol cause marked flushing and nausea."},
            {"when": {"rsid": "rs671", "genotype": ["AG"]},
             "call": "Alcohol flush reaction likely",
             "detail": "One ALDH2*2 copy reduces acetaldehyde clearance sharply, which is the classic flushing response."},
            {"when": {"rsid": "rs671", "genotype": ["GG"]},
             "call": "No flush reaction expected",
             "detail": "Both ALDH2 copies are fully functional, so acetaldehyde is cleared at the usual rate."},
        ],
        "default_call": "Undetermined",
        "magnitude": 2.5,
        "evidence": "ALDH2 Glu504Lys is one of the largest single-SNP effects known for a common trait, and is near-monomorphic outside East Asia.",
        "caveat": "The variant allele is common in East Asian ancestry and rare elsewhere, so most non-East-Asian results will simply be GG. Flushing severity also varies with drinking pattern and ADH1B genotype, which this trait does not include.",
    },
    {
        "key": "caffeine_metabolism",
        "name": "Caffeine metabolism speed",
        "category": "Metabolism",
        "rsids": ["rs762551"],
        "rules": [
            {"when": {"rsid": "rs762551", "genotype": ["AA"]},
             "call": "Fast caffeine metabolism",
             "detail": "CYP1A2*1F homozygote, associated with quicker caffeine clearance and shorter perceived effect."},
            {"when": {"rsid": "rs762551", "genotype": ["AC"]},
             "call": "Intermediate caffeine metabolism",
             "detail": "One fast and one slow allele, giving an intermediate clearance rate."},
            {"when": {"rsid": "rs762551", "genotype": ["CC"]},
             "call": "Slow caffeine metabolism",
             "detail": "Slower CYP1A2 induction, so caffeine and several CYP1A2 substrates linger longer."},
        ],
        "default_call": "Undetermined",
        "magnitude": 1.5,
        "evidence": "CYP1A2 rs762551 is the standard caffeine clearance marker and appears in most caffeine pharmacokinetic studies.",
        "caveat": "Smoking, oral contraceptives and heavy coffee drinking all induce or inhibit CYP1A2 more strongly than this genotype does, so measured clearance often disagrees with the predicted class.",
    },
    {
        "key": "bitter_taste_ptc",
        "name": "Bitter taste perception (PTC)",
        "category": "Taste",
        "rsids": ["rs713598", "rs1726866"],
        "rules": [
            {"when": {"all": [{"rsid": "rs713598", "genotype": ["GG"]},
                              {"rsid": "rs1726866", "min_copies": {"allele": "A", "n": 1}}]},
             "call": "Strong bitter taster",
             "detail": "Two taster copies at the main TAS2R38 site with a taster allele at the second site, the classic PAV/PAV super-taster pattern."},
            {"when": {"rsid": "rs713598", "min_copies": {"allele": "G", "n": 1}},
             "call": "Bitter taster",
             "detail": "At least one TAS2R38 taster allele, so PTC and related compounds in brassicas usually taste bitter."},
            {"when": {"rsid": "rs713598", "genotype": ["CC"]},
             "call": "Bitter non-taster",
             "detail": "Two non-taster copies, so PTC-type bitterness is much weaker or absent."},
        ],
        "default_call": "Undetermined",
        "magnitude": 1.5,
        "evidence": "TAS2R38 haplotypes built from rs713598, rs1726866 and rs10246939 explain most of the variance in PTC tasting.",
        "caveat": "The full effect depends on the three-site haplotype and this trait reads only two of the three sites, so intermediate tasters are not separated cleanly. Perceived bitterness of real food also depends on the compound, not just PTC.",
    },
    {
        "key": "asparagus_odour_detection",
        "name": "Asparagus urine odour detection",
        "category": "Smell",
        "rsids": ["rs4481887"],
        "rules": [
            {"when": {"rsid": "rs4481887", "min_copies": {"allele": "G", "n": 1}},
             "call": "Can probably detect asparagus urine odour",
             "detail": "Carries the allele associated with being able to smell the sulphur metabolites of asparagus."},
            {"when": {"rsid": "rs4481887", "genotype": ["AA"]},
             "call": "Reduced ability to detect asparagus urine odour",
             "detail": "Associated with specific anosmia for the asparagus metabolites, so the odour may not register at all."},
        ],
        "default_call": "Undetermined",
        "magnitude": 1.0,
        "evidence": "A large replicated olfactory receptor cluster association on chromosome 1 for asparagus anosmia.",
        "caveat": "This is a single tag in a dense olfactory receptor cluster, effect sizes are modest, and everyone produces the metabolites regardless of whether they can smell them.",
    },
]

TRAITS += [
    {
        "key": "cilantro_preference",
        "name": "Cilantro (coriander) preference",
        "category": "Taste",
        "rsids": ["rs72921001"],
        "rules": [
            {"when": {"rsid": "rs72921001", "min_copies": {"allele": "C", "n": 1}},
             "call": "More likely to find cilantro soapy",
             "detail": "Carries the OR6A2-linked allele associated with perceiving the aldehydes in cilantro as soapy."},
            {"when": {"rsid": "rs72921001", "genotype": ["AA"]},
             "call": "No soapy-taste association",
             "detail": "Does not carry the soapy-taste allele, so dislike of cilantro is more likely to be learned than genetic."},
        ],
        "default_call": "Undetermined",
        "magnitude": 1.0,
        "evidence": "Replicated association near the OR6A2 olfactory receptor with self-reported soapy cilantro taste.",
        "caveat": "The genetic effect on cilantro liking is small and explains only a few percent of the variance. Culture and exposure matter far more, so plenty of AA people dislike cilantro anyway.",
    },
    {
        "key": "earwax_and_body_odour",
        "name": "Earwax type and underarm odour",
        "category": "Physical",
        "rsids": ["rs17822931"],
        "rules": [
            {"when": {"rsid": "rs17822931", "genotype": ["TT"]},
             "call": "Dry earwax and reduced underarm odour",
             "detail": "Two non-functional ABCC11 copies give dry, flaky earwax and much less axillary odour precursor secretion."},
            {"when": {"rsid": "rs17822931", "genotype": ["CT"]},
             "call": "Wet earwax, carrier of the dry allele",
             "detail": "One functional ABCC11 copy is enough for wet earwax and typical odour."},
            {"when": {"rsid": "rs17822931", "genotype": ["CC"]},
             "call": "Wet earwax and typical underarm odour",
             "detail": "Two functional ABCC11 copies, the common pattern in European and African ancestry."},
        ],
        "default_call": "Undetermined",
        "magnitude": 2.0,
        "evidence": "ABCC11 Gly180Arg is effectively deterministic for earwax type, one of the cleanest genotype to phenotype maps in humans.",
        "caveat": "Earwax type is close to determined by this variant, but body odour is only partly explained by it: skin microbiome, diet and hygiene all contribute, so a dry-earwax result does not guarantee no odour.",
    },
    {
        "key": "muscle_fibre_actn3",
        "name": "Muscle fibre type (ACTN3 R577X)",
        "category": "Fitness",
        "rsids": ["rs1815739"],
        "rules": [
            {"when": {"rsid": "rs1815739", "genotype": ["CC"]},
             "call": "Power-leaning (two functional ACTN3 copies)",
             "detail": "Two R577 copies produce alpha-actinin-3 in fast-twitch fibres, over-represented among elite sprint and power athletes."},
            {"when": {"rsid": "rs1815739", "genotype": ["CT"]},
             "call": "Mixed power and endurance profile",
             "detail": "One functional copy, the most common genotype, with no strong lean either way."},
            {"when": {"rsid": "rs1815739", "genotype": ["TT"]},
             "call": "Endurance-leaning (no alpha-actinin-3)",
             "detail": "Two stop-codon copies mean no alpha-actinin-3 at all, which is common in the general population and over-represented in endurance athletes."},
        ],
        "default_call": "Undetermined",
        "magnitude": 1.5,
        "evidence": "ACTN3 R577X is the most studied fitness variant, with consistent athlete-cohort differences.",
        "caveat": "Effects on ordinary performance are small and training beats genotype comfortably. Roughly a fifth of people worldwide are XX including successful sprinters, so this is a lean, not a limit.",
    },
    {
        "key": "photic_sneeze_reflex",
        "name": "Photic sneeze reflex",
        "category": "Reflex",
        "rsids": ["rs10427255"],
        "rules": [
            {"when": {"rsid": "rs10427255", "min_copies": {"allele": "C", "n": 1}},
             "call": "More likely to sneeze in bright light",
             "detail": "Carries the allele associated with the photic sneeze reflex, sneezing on sudden exposure to bright sunlight."},
            {"when": {"rsid": "rs10427255", "genotype": ["TT"]},
             "call": "Less likely to sneeze in bright light",
             "detail": "Does not carry the associated allele."},
        ],
        "default_call": "Undetermined",
        "magnitude": 1.0,
        "evidence": "Discovered and replicated in large self-reported cohort studies of the photic sneeze reflex.",
        "caveat": "Based on self-reported phenotypes with a modest odds ratio, so it shifts probability rather than predicting the reflex. Many carriers never notice it.",
    },
    {
        "key": "mc1r_freckling",
        "name": "Freckling and red hair (MC1R)",
        "category": "Appearance",
        "rsids": ["rs1805007", "rs1805008"],
        "rules": [
            {"when": {"any": [{"rsid": "rs1805007", "genotype": ["TT"]},
                              {"rsid": "rs1805008", "genotype": ["TT"]}]},
             "call": "Two MC1R red-hair variants",
             "detail": "Homozygous for a strong MC1R loss-of-function allele, which makes red hair, heavy freckling and poor tanning very likely."},
            {"when": {"all": [{"rsid": "rs1805007", "min_copies": {"allele": "T", "n": 1}},
                              {"rsid": "rs1805008", "min_copies": {"allele": "T", "n": 1}}]},
             "call": "Compound MC1R variant carrier",
             "detail": "One copy of each of the two strong MC1R alleles, a combination that behaves much like a red-hair genotype."},
            {"when": {"any": [{"rsid": "rs1805007", "min_copies": {"allele": "T", "n": 1}},
                              {"rsid": "rs1805008", "min_copies": {"allele": "T", "n": 1}}]},
             "call": "MC1R variant carrier",
             "detail": "One strong MC1R allele, associated with freckling, sun sensitivity and a raised chance of red hair in children."},
        ],
        "default_call": "No MC1R R151C or R160W variant detected",
        "magnitude": 2.0,
        "evidence": "MC1R R151C (rs1805007) and R160W (rs1805008) are the two strongest common red-hair alleles.",
        "caveat": "MC1R has many other loss-of-function alleles, including R142H and D294H, that these two SNPs do not cover, so a negative result does not rule out red hair or high sun sensitivity. Sun protection advice does not depend on genotype.",
    },
    {
        "key": "eye_colour_herc2",
        "name": "Eye colour (HERC2/OCA2)",
        "category": "Appearance",
        "rsids": ["rs12913832"],
        "rules": [
            {"when": {"rsid": "rs12913832", "genotype": ["GG"]},
             "call": "Blue or light eyes likely",
             "detail": "Two copies of the low-OCA2-expression allele, the pattern found in the large majority of blue-eyed Europeans."},
            {"when": {"rsid": "rs12913832", "genotype": ["AG"]},
             "call": "Intermediate or brown eyes",
             "detail": "One copy each way, which gives the widest spread of outcomes: green, hazel and brown are all common."},
            {"when": {"rsid": "rs12913832", "genotype": ["AA"]},
             "call": "Brown eyes likely",
             "detail": "Two copies of the high-expression allele, which normally gives brown eyes."},
        ],
        "default_call": "Undetermined",
        "magnitude": 2.0,
        "evidence": "HERC2 rs12913832 is the single largest-effect common eye colour variant, acting on OCA2 expression.",
        "caveat": "This variant explains a large share of blue versus brown eye colour in Europeans and very little outside Europe. It is a probability, not a determination: intermediate genotypes especially can produce almost any eye colour, and many genes contribute.",
    },
]

TRAITS += [
    {
        "key": "milk_digestion_summary",
        "name": "Lactose and milk digestion summary",
        "category": "Diet",
        "rsids": ["rs4988235", "rs182549"],
        "rules": [
            {"when": {"all": [{"rsid": "rs4988235", "genotype": ["CC"]},
                              {"rsid": "rs182549", "genotype": ["CC"]}]},
             "call": "Both LCT tags non-persistent",
             "detail": "Neither MCM6 tag carries a persistence allele, so lactase activity most likely declined after childhood and fresh milk is often poorly tolerated. Aged cheese, yoghurt and lactase supplements usually still work."},
            {"when": {"all": [{"rsid": "rs4988235", "min_copies": {"allele": "T", "n": 1}},
                              {"rsid": "rs182549", "min_copies": {"allele": "T", "n": 1}}]},
             "call": "Both LCT tags show persistence",
             "detail": "Both tags carry a persistence allele, the usual European lactase persistence haplotype, so milk is normally digested throughout life."},
            {"when": {"any": [{"rsid": "rs4988235", "min_copies": {"allele": "T", "n": 1}},
                              {"rsid": "rs182549", "min_copies": {"allele": "T", "n": 1}}]},
             "call": "At least one persistence allele",
             "detail": "One tag carries a persistence allele, which is usually enough for adult lactose digestion even if tolerance is not unlimited."},
            {"when": {"rsid": "rs4988235", "genotype": ["CC"]},
             "call": "Primary LCT tag non-persistent",
             "detail": "The main tag shows no persistence allele; the second tag was not available to confirm the haplotype."},
        ],
        "default_call": "Undetermined",
        "magnitude": 1.5,
        "evidence": "rs4988235 and rs182549 sit on the same MCM6 enhancer haplotype and are near-perfectly correlated in Europeans.",
        "caveat": "Both tags describe the European haplotype only, and neither predicts symptom severity. Self-reported dairy intolerance frequently has nothing to do with lactase, so a breath test or an elimination trial is more informative than genotype.",
    },
    {
        "key": "nicotine_dependence_liability",
        "name": "Nicotine dependence liability",
        "category": "Behaviour",
        "rsids": ["rs1051730"],
        "rules": [
            {"when": {"rsid": "rs1051730", "genotype": ["AA"]},
             "call": "Higher nicotine dependence liability",
             "detail": "Two copies of the CHRNA3 risk allele, associated with smoking more cigarettes per day and finding cessation harder among people who do smoke."},
            {"when": {"rsid": "rs1051730", "genotype": ["AG"]},
             "call": "Moderately higher nicotine dependence liability",
             "detail": "One copy of the risk allele, with an intermediate effect on cigarettes per day."},
            {"when": {"rsid": "rs1051730", "genotype": ["GG"]},
             "call": "Baseline nicotine dependence liability",
             "detail": "No copy of the risk allele at this position."},
        ],
        "default_call": "Undetermined",
        "magnitude": 1.5,
        "evidence": "CHRNA3 rs1051730 is one of the most replicated smoking quantity loci, with a consistent effect of roughly one extra cigarette per day per allele.",
        "caveat": "This predicts intensity among smokers, not whether anyone starts smoking, and it has no meaning for a person who has never smoked. The effect is about one cigarette per day per allele, which is tiny next to social and environmental drivers.",
    },
    {
        "key": "chronotype_clock",
        "name": "Morning versus evening chronotype",
        "category": "Sleep",
        "rsids": ["rs1801260"],
        "rules": [
            {"when": {"rsid": "rs1801260", "genotype": ["GG"]},
             "call": "Evening-leaning chronotype",
             "detail": "Two copies of the CLOCK 3111 variant allele, associated on average with later bedtimes and later peak alertness."},
            {"when": {"rsid": "rs1801260", "genotype": ["AG"]},
             "call": "Intermediate chronotype",
             "detail": "One copy of the variant allele, with no clear lean."},
            {"when": {"rsid": "rs1801260", "genotype": ["AA"]},
             "call": "Morning-leaning chronotype",
             "detail": "No copy of the variant allele, associated on average with earlier sleep timing."},
        ],
        "default_call": "Undetermined",
        "magnitude": 1.0,
        "evidence": "CLOCK 3111T/C is the classic candidate chronotype variant, though modern GWAS place most chronotype variance elsewhere.",
        "caveat": "Vendors report this position on different strands, so the same person can appear to have opposite genotypes in two files. The effect is also small: large chronotype GWAS find hundreds of loci, and shift work, light exposure and age move sleep timing far more.",
    },
    {
        "key": "ace_endurance_power",
        "name": "Endurance versus power (ACE)",
        "category": "Fitness",
        "rsids": ["rs4343"],
        "rules": [
            {"when": {"rsid": "rs4343", "genotype": ["AA"]},
             "call": "Insertion-like, endurance-leaning",
             "detail": "Proxy genotype for the ACE I/I insertion pattern, over-represented in endurance athlete cohorts."},
            {"when": {"rsid": "rs4343", "genotype": ["AG"]},
             "call": "Mixed insertion and deletion pattern",
             "detail": "Proxy genotype for ACE I/D, the most common pattern, with no strong lean."},
            {"when": {"rsid": "rs4343", "genotype": ["GG"]},
             "call": "Deletion-like, power-leaning",
             "detail": "Proxy genotype for ACE D/D, associated with higher ACE activity and with power and sprint phenotypes."},
        ],
        "default_call": "Undetermined",
        "magnitude": 1.0,
        "evidence": "rs4343 is the standard SNP proxy for the ACE intron 16 insertion/deletion polymorphism, which arrays cannot type directly.",
        "caveat": "This is a proxy for an insertion/deletion the array does not genotype, so the mapping is imperfect and ancestry dependent. Athlete-cohort associations for ACE have replicated inconsistently and say almost nothing about an individual.",
    },
    {
        "key": "pain_sensitivity_faah",
        "name": "Pain sensitivity and anxiety (FAAH)",
        "category": "Neurology",
        "rsids": ["rs324420"],
        "rules": [
            {"when": {"rsid": "rs324420", "genotype": ["AA"]},
             "call": "Reduced FAAH activity",
             "detail": "Two copies of the C385A variant reduce fatty acid amide hydrolase activity, raising anandamide tone, which is associated with lower reported anxiety and pain sensitivity."},
            {"when": {"rsid": "rs324420", "genotype": ["AC"]},
             "call": "Intermediate FAAH activity",
             "detail": "One variant copy, with an intermediate effect on endocannabinoid tone."},
            {"when": {"rsid": "rs324420", "genotype": ["CC"]},
             "call": "Typical FAAH activity",
             "detail": "No variant copy, so endocannabinoid breakdown proceeds at the usual rate."},
        ],
        "default_call": "Undetermined",
        "magnitude": 1.0,
        "evidence": "FAAH C385A reliably reduces enzyme stability in vitro, with consistent effects on anandamide levels.",
        "caveat": "The enzyme effect is solid but the behavioural associations are small and inconsistent across studies. Pain sensitivity is highly context dependent and this genotype must not be used to judge anyone's reported pain.",
    },
    {
        "key": "oxytocin_social_sensitivity",
        "name": "Oxytocin receptor social sensitivity",
        "category": "Social",
        "rsids": ["rs53576"],
        "rules": [
            {"when": {"rsid": "rs53576", "genotype": ["GG"]},
             "call": "Higher reported social sensitivity",
             "detail": "GG homozygotes report slightly higher empathy and prosocial behaviour on average in several studies."},
            {"when": {"rsid": "rs53576", "genotype": ["AG"]},
             "call": "Intermediate reported social sensitivity",
             "detail": "One A allele, with an intermediate average score."},
            {"when": {"rsid": "rs53576", "genotype": ["AA"]},
             "call": "Lower reported social sensitivity",
             "detail": "AA homozygotes score slightly lower on average on empathy measures in the same studies."},
        ],
        "default_call": "Undetermined",
        "magnitude": 0.5,
        "evidence": "OXTR rs53576 is the most studied social behaviour candidate SNP, though the literature is mixed.",
        "caveat": "This is a small candidate-gene effect from an area with well-documented replication problems, and effects are strongly context dependent. It says nothing about any individual's empathy and must never be read as a personality diagnosis.",
    },
    {
        "key": "fat_taste_cd36",
        "name": "Fat taste perception (CD36)",
        "category": "Taste",
        "rsids": ["rs1761667"],
        "rules": [
            {"when": {"rsid": "rs1761667", "genotype": ["AA"]},
             "call": "Reduced oral fat sensitivity",
             "detail": "Associated with lower CD36 expression and a higher detection threshold for fatty acids, which tends to go with a preference for fattier food."},
            {"when": {"rsid": "rs1761667", "genotype": ["AG"]},
             "call": "Intermediate oral fat sensitivity",
             "detail": "One copy each way, with an intermediate detection threshold."},
            {"when": {"rsid": "rs1761667", "genotype": ["GG"]},
             "call": "Higher oral fat sensitivity",
             "detail": "Associated with higher CD36 expression and a lower fatty acid detection threshold."},
        ],
        "default_call": "Undetermined",
        "magnitude": 0.5,
        "evidence": "CD36 rs1761667 shows reproducible differences in oral fatty acid detection thresholds in sensory studies.",
        "caveat": "Sensory studies use small panels and the effect on real-world food choice is weak. Habitual diet changes fat taste thresholds within weeks, which swamps the genotype effect.",
    },
]

TRAIT_KEYS = tuple(t["key"] for t in TRAITS)


# ---------------------------------------------------------------------------
# Rule evaluator
# ---------------------------------------------------------------------------

def _spec_rsids(spec: Any) -> list[str]:
    """Collect every rsID mentioned by a rule spec, in declaration order."""
    found: list[str] = []
    if not isinstance(spec, dict):
        return found
    for combinator in ("all", "any"):
        for sub in (spec.get(combinator) or []):
            for rsid in _spec_rsids(sub):
                if rsid not in found:
                    found.append(rsid)
    rsid = str(spec.get("rsid") or "").strip().lower()
    if rsid and rsid not in found:
        found.append(rsid)
    return found


def _eval_spec(spec: Any, genotypes: dict) -> bool | None:
    """Evaluate one rule spec against a genotype mapping.

    Returns True or False when the spec could be decided, and None when a
    required rsID was not genotyped, so the caller can fall through to the
    next rule instead of treating "unknown" as "false".
    """
    if not isinstance(spec, dict):
        return None
    if "all" in spec:
        results = [_eval_spec(sub, genotypes) for sub in (spec.get("all") or [])]
        if not results or any(r is None for r in results):
            return None
        return all(results)
    if "any" in spec:
        results = [_eval_spec(sub, genotypes) for sub in (spec.get("any") or [])]
        if not results:
            return None
        if any(r is True for r in results):
            return True
        if any(r is None for r in results):
            return None
        return False
    rsid = str(spec.get("rsid") or "").strip().lower()
    if not rsid:
        return None
    key = _read(genotypes, rsid)
    if key is None or key == NOCALL_KEY:
        return None
    if "genotype" in spec:
        wanted = {_key_from_text(text) for text in (spec.get("genotype") or [])}
        return key in wanted
    if "min_copies" in spec:
        rule = spec.get("min_copies") or {}
        allele = str(rule.get("allele") or "").upper()
        needed = int(rule.get("n", 1))
        return sum(1 for base in key if base == allele) >= needed
    if "max_copies" in spec:
        rule = spec.get("max_copies") or {}
        allele = str(rule.get("allele") or "").upper()
        allowed = int(rule.get("n", 0))
        return sum(1 for base in key if base == allele) <= allowed
    return True


def evaluate_trait(trait: dict, genotypes: dict) -> dict | None:
    """Evaluate one trait definition against a genotype mapping.

    Returns None when none of the trait's rsIDs produced a real call, so the
    report never shows a trait the file cannot speak to. Otherwise the first
    rule that evaluates True supplies the call and detail; if no rule matches,
    the trait's default_call is used.
    """
    g = _index(genotypes)
    rsids = [str(r).strip().lower() for r in (trait.get("rsids") or [])]
    called = [rsid for rsid in rsids if _has_call(g, rsid)]
    if not called:
        return None

    call = trait.get("default_call", "Undetermined")
    detail = ""
    decided = called[0]
    for rule in (trait.get("rules") or []):
        if _eval_spec(rule.get("when"), g) is True:
            call = rule.get("call", call)
            detail = rule.get("detail", "")
            picked = [r for r in _spec_rsids(rule.get("when")) if _has_call(g, r)]
            if picked:
                decided = picked[0]
            break

    coverage = round(len(called) / len(rsids), 3) if rsids else 0.0
    return {
        "key": trait.get("key", ""),
        "name": trait.get("name", ""),
        "category": trait.get("category", ""),
        "call": call,
        "detail": detail,
        "genotype": _read(g, decided) or NOCALL_KEY,
        "rsid": decided,
        "magnitude": float(trait.get("magnitude", 0.0)),
        "evidence": trait.get("evidence", ""),
        "caveat": trait.get("caveat", ""),
        "coverage": coverage,
    }


def predict_traits(genotypes: dict) -> list[dict]:
    """Evaluate every trait in TRAITS, skipping those the file cannot answer."""
    results: list[dict] = []
    for trait in TRAITS:
        result = evaluate_trait(trait, genotypes)
        if result is not None:
            results.append(result)
    return results


def _zygosity(key: str) -> str:
    """Classify a genotype key as homozygous, heterozygous or no_call."""
    if not key or key == NOCALL_KEY or len(key) != 2:
        return "no_call"
    return "homozygous" if key[0] == key[1] else "heterozygous"


def to_findings(trait_results: list[dict], blood: dict | None = None) -> list[dict]:
    """Convert trait results (and optionally a blood type result) to findings.

    Traits are neutral. Every finding produced here has entity_type "trait",
    silo "informational" and an empty repute, so the UI can never colour a
    trait as good or bad. The finding's rsid is the trait key rather than a
    single rsID, because several traits read more than one position.
    """
    findings: list[dict] = []
    for result in (trait_results or []):
        genotype = result.get("genotype", NOCALL_KEY)
        allele1 = genotype[0] if len(genotype) == 2 else "N"
        allele2 = genotype[1] if len(genotype) == 2 else "N"
        detail = result.get("detail", "")
        caveat = result.get("caveat", "")
        summary = result.get("call", "")
        if detail:
            summary = f"{summary}. {detail}"
        findings.append({
            "rsid": result.get("key", ""),
            "entity_type": "trait",
            "name": result.get("name", ""),
            "gene": "",
            "chromosome": "",
            "position": 0,
            "allele1": allele1,
            "allele2": allele2,
            "genotype": genotype,
            "zygosity": _zygosity(genotype),
            "clinical_sig": "informational",
            "conditions": result.get("name", ""),
            "summary": summary,
            "interpretation": " ".join(part for part in (
                detail, result.get("evidence", ""), f"Limitation: {caveat}" if caveat else ""
            ) if part),
            "category": result.get("category", ""),
            "silo": "informational",
            "repute": "",
            "magnitude": float(result.get("magnitude", 0.0)),
            "sources": ["dnainsight_traits"],
            "caveat": caveat,
            "coverage": result.get("coverage", 0.0),
            "source_rsid": result.get("rsid", ""),
        })

    if blood:
        abo = blood.get("abo") or {}
        rh = blood.get("rh") or {}
        findings.append({
            "rsid": "blood_type",
            "entity_type": "trait",
            "name": "ABO and RhD blood group",
            "gene": "ABO/RHD",
            "chromosome": "",
            "position": 0,
            "allele1": "N",
            "allele2": "N",
            "genotype": abo.get("genotype_call", "") or NOCALL_KEY,
            "zygosity": "no_call",
            "clinical_sig": "informational",
            "conditions": "Blood group",
            "summary": blood.get("summary", ""),
            "interpretation": " ".join(part for part in (
                blood.get("summary", ""),
                f"ABO limitation: {abo.get('caveat', '')}" if abo.get("caveat") else "",
                f"RhD limitation: {rh.get('caveat', '')}" if rh.get("caveat") else "",
            ) if part),
            "category": "Blood",
            "silo": "informational",
            "repute": "",
            "magnitude": 1.0,
            "sources": ["dnainsight_traits"],
            "caveat": abo.get("caveat", ""),
            "coverage": 1.0 if blood.get("blood_type") != "unknown" else 0.0,
            "source_rsid": ABO_DECISIVE_TAG,
        })
    return findings
