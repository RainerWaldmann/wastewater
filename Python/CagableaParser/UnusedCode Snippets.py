"""VariantsNew"""

for voc in vg.variantdict:
    if voc not in detailedCounts:  # if no data generate 0 mean data
        for sample in samples:
            variantDataMeansSD[sample].loc[voc] = VariantData(0, np.nan, 0, "", True)
    else:
        for sample in samples:
            x = detailedCounts[voc][sample]  # Var freqs for one variant for one sample
            mask = [a == a for a in x]  # nans are marked false
            hasRequiredMuts = True
            if vg.variantdict[
                voc].minstarredmutsforpass > 0:  # and len(vg.variantdict[voc].getRequiredMutations()) != 0 :
                # requiredMuts = [b for b in x.index if any([a in b for a in vg.variantdict[voc].getRequiredMutations()])]
                requiredMutsFound = x[mask][
                    [b for b in x[mask].index if any([a in b for a in vg.variantdict[voc].getRequiredMutations()])]]
                hasRequiredMuts = len(requiredMutsFound) >= vg.variantdict[voc].minstarredmutsforpass and all(
                    [c > 0 for c in requiredMutsFound])

            # mask = [True] * len(x)
            # nans = [a != a for a in x] #nans are marked true

            # m = x.mean(skipna=True)
            # sd = np.nan if x.count() < 2 else x.std(skipna=True)
            if hasRequiredMuts:
                if sum(mask) >= 3:  # mask outliers
                    mask = np.logical_and(getOutlierMask(x[mask]), mask)
                    #
                    # #outlier_mask = x.between(m - fold_sd_dev_from_mean * sd,
                    # #                        m + fold_sd_dev_from_mean * sd)  # mask with outliers
                    # outlier_mask = getOutlierMask(x)
                    # doneMuts_mask = np.logical_or(nans,outlier_mask)  # will set to true if either NaN or not an outlier
                    # while (sum(mask) / len(mask)) >= 0.60 and (sum(mask) >= 2) \
                    #         and (sum(~doneMuts_mask) > 0): # positions not finished must remain # at least 70% of muts and at least two muts not maked
                    #     diff_from_mean = [abs(v - m) for v in x]  # list wit differences between value and mean
                    #     # find index of biggest outlier
                    #     biggest_outlier_index = []  # holds list of indices of biggest outliers
                    #     biggest_outlier_value = -1
                    #     for i in range(len(x)):
                    #         if not outlier_mask[i] and not doneMuts_mask[i]: # is an outlier and not done yet
                    #             if diff_from_mean[i] > biggest_outlier_value:
                    #                 biggest_outlier_value = diff_from_mean[i]
                    #                 biggest_outlier_index = [i]
                    #             elif diff_from_mean[i] == biggest_outlier_value:  # outlier has same dev as previous
                    #                 biggest_outlier_index.append(i)
                    #     for i in range(len(biggest_outlier_index)):
                    #         mask[biggest_outlier_index[i]] = False
                    #         doneMuts_mask[biggest_outlier_index[i]] = True
                    #     m = x[mask].mean(skipna=True)
                    #     sd = x[mask].std(skipna=True)
                    #     outlier_mask = getOutlierMask(x[mask])
                    #     #outlier_mask = x.between(m - fold_sd_dev_from_mean * sd, m + fold_sd_dev_from_mean * sd)
                    #     doneMuts_mask = np.logical_or(doneMuts_mask, outlier_mask)  # set all that pass to True

            if sum(mask) > 1:
                box = getBoxPlot(x, x[mask]).to_html(full_html=False, include_plotlyjs='cdn') if \
                    settings.plotBoxPlotsForDetailedVarCounts else ""
            else:
                box = ""

            x = x[mask]
            detailedCountsMask[voc][sample] = mask
            if sum(mask) > 0 and hasRequiredMuts:  # at least one left
                variantDataMeansSD[sample].loc[voc] = VariantData(x.mean(skipna=True),
                                                                  np.nan if sum(mask) <= 1 else x.std(
                                                                      skipna=True) / math.sqrt(sum(mask)),
                                                                  len(x), box, False)
            else:
                variantDataMeansSD[sample].loc[voc] = VariantData(0, np.nan, 0, "", True)
"""xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"""
# # populate variantDataMeansSD
# for voc in vg.variantdict:
#     if voc not in detailedCounts:  # if no data generate 0 mean data
#         for sample in samples:
#             variantDataMeansSD[sample].loc[voc] = VariantData(0, np.nan, 0, "", True)
#     else:
#         for sample in samples:
#             x = detailedCounts[voc][sample]
#             mask = [True] * len(x)  # fill mask (info about outliers) with TRUE
#             m = x.mean(skipna=True)
#             sd = np.nan if x.count() < 2 else x.std(skipna=True)
#             if len(x) >= 3:  # mask outliers
#                 nans = [a != a for a in x]
#                 outlier_mask = x.between(m - fold_sd_dev_from_mean * sd,
#                                         m + fold_sd_dev_from_mean * sd)  # mask with outliers
#                 msk_final = [True] * len(x)
#                 doneMuts_mask = np.logical_or(outlier_mask,
#                                               nans)  # will set to true if done for each false (masked) or NAN pos
#                 while (sum(msk_final) / len(msk_final)) >= 0.60 and (sum(msk_final) >= 2) and (
#                         sum(~doneMuts_mask) > 0):
#                     diff_from_mean = [abs(v - m) for v in x]  # difference between value and mean
#                     # find index of biggest outlier for considered - dont_consider = FALSE
#                     biggest_outlier_index = []  # holds list of indices of biggest outliers
#                     biggest_outlier_value = 0
#                     for i in range(len(x)):
#                         if not outlier_mask[i] and not doneMuts_mask[i]:
#                             if diff_from_mean[i] > biggest_outlier_value:
#                                 biggest_outlier_value = diff_from_mean[i]
#                                 biggest_outlier_index = [i]
#                             elif diff_from_mean[i] == biggest_outlier_value:  # outlier has same dev as previous
#                                 biggest_outlier_index.append(i)
#                     for i in range(len(biggest_outlier_index)):
#                         msk_final[biggest_outlier_index[i]] = False
#                         doneMuts_mask[biggest_outlier_index[i]] = True
#                     m = x[msk_final].mean(skipna=True)
#                     sd = x[msk_final].std(skipna=True)
#                     outlier_mask = x.between(m - fold_sd_dev_from_mean * sd, m + fold_sd_dev_from_mean * sd)
#                     doneMuts_mask = np.logical_or(doneMuts_mask, outlier_mask)  # set all that pass to True
#
#                 mask = np.logical_and(msk_final, np.invert(nans))
#
#             if x.count() > 1:
#                 box = getBoxPlot(x, x[mask]).to_html(full_html=False, include_plotlyjs='cdn') if \
#                     settings.plotBoxPlotsForDetailedVarCounts else ""
#             else:
#                 box = ""
#
#             x = x[mask]
#             detailedCountsMask[voc][sample] = mask
#             variantDataMeansSD[sample].loc[voc] = VariantData(x.mean(skipna=True),
#                                                               np.nan if len(x) <= 1 else x.std(
#                                                                   skipna=True) / math.sqrt(len(x)),
#                                                               len(x), box, False)

"""xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"""