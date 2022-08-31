import numpy as np

def Binning_distance_from_ice_edge(ice_achv, var_achv, ntime, nlat, nlon, hem, outputcount = False, model = False):
    for j in range(ntime):
        if model:
            ice = ice_achv[j, :, :]
        else:
            ice = ice_achv[:, :, j]
        var = var_achv[j]
        ndim = var.ndim

        for i in range(nlon):
            tmp = ice[:,i]
            if ndim==3:
                var_tmp = var[:,:,i]
                nh = var_tmp.shape[0]
            else:
                var_tmp = var[:,i]
            
            if hem == 'nh':
                ind1 = np.where(tmp <= 0.5)[0]
                ind2 = np.where(tmp > 0.5)[0]
            else:
                ind1 = np.where(tmp > 0.5)[0]
                ind2 = np.where(tmp <= 0.5)[0]                

            n1 = ind1.size
            n2 = ind2.size

            if i == 0:
                tmp0 = np.full((2*nlat), -999.)
                tmp0[nlat - n1: nlat] = tmp[ind1]
                tmp0[nlat: nlat + n2] = tmp[ind2]

                count0 = np.full((2*nlat), 1.)
                
                if ndim==3:
                    var_tmp0 = np.full((nh, 2*nlat), -999.)
                    var_tmp0[:,nlat - n1: nlat] = var_tmp[:,ind1]
                    var_tmp0[:,nlat: nlat + n2] = var_tmp[:,ind2]
                else:
                    var_tmp0 = np.full((2*nlat), -999.)
                    var_tmp0[nlat - n1: nlat] = var_tmp[ind1]
                    var_tmp0[nlat: nlat + n2] = var_tmp[ind2]
            else:
                tmp1 = np.full((2*nlat), -999.)
                tmp1[nlat - n1: nlat] = tmp[ind1]
                tmp1[nlat: nlat + n2] = tmp[ind2]

                count1 = np.full((2*nlat), 1.)
                
                if ndim==3:
                    var_tmp1 = np.full((nh, 2*nlat), -999.)
                    var_tmp1[:,nlat - n1: nlat] = var_tmp[:,ind1]
                    var_tmp1[:,nlat: nlat + n2] = var_tmp[:,ind2]
                else:
                    var_tmp1 = np.full((2*nlat), -999.)
                    var_tmp1[nlat - n1: nlat] = var_tmp[ind1]
                    var_tmp1[nlat: nlat + n2] = var_tmp[ind2]

                tmp0 = np.dstack((tmp0, tmp1))
                count0 = np.dstack((count0, count1))
                var_tmp0 = np.dstack((var_tmp0, var_tmp1))

        var_tmp0[var_tmp0 == -999.] = np.nan
        tmp0[tmp0 == -999.] = np.nan
        
        if ndim == 2:
            count0[np.isnan(var_tmp0)] = 0.
        
        count00 = nlon*np.nanmean(count0[0], axis = 1)
        tmp00 = np.nanmean(tmp0[0], axis = 1)
        if ndim == 3:
            var_tmp00 = np.nanmean(var_tmp0, axis = 2)
        else:
            var_tmp00 = np.nanmean(var_tmp0[0], axis = 1)

        if j == 0:
            count_out = count00
            var_out = var_tmp00
            ice_out = tmp00
        else:
            count_out = np.dstack((count_out, count00))        
            var_out = np.dstack((var_out, var_tmp00))
            ice_out = np.dstack((ice_out, tmp00))

    count_out = count_out[0]
    ice_out = ice_out[0]
    if ndim != 3:
        var_out = var_out[0]
    
    if ndim == 3:
        return var_out, ice_out
    else:
        if outputcount:
            return var_out, ice_out, count_out    
        else:
            return var_out, ice_out

# def Binning_distance_from_ice_edge_model(ice_achv, var_achv, ntime, nlat, nlon, hem):
#     for j in range(ntime):
#         ice = ice_achv[j, :, :]
#         var = var_achv[j]
#         ndim = var.ndim

#         for i in range(nlon):
#             tmp = ice[:,i]
#             if ndim==3:
#                 var_tmp = var[:,:,i]
#                 nh = var_tmp.shape[0]
#             else:
#                 var_tmp = var[:,i]
            
#             if hem == 'nh':
#                 ind1 = np.where(tmp <= 0.5)[0]
#                 ind2 = np.where(tmp > 0.5)[0]
#             else:
#                 ind1 = np.where(tmp > 0.5)[0]
#                 ind2 = np.where(tmp <= 0.5)[0]                

#             n1 = ind1.size
#             n2 = ind2.size

#             if i == 0:
#                 tmp0 = np.full((2*nlat), -999.)
#                 tmp0[nlat - n1: nlat] = tmp[ind1]
#                 tmp0[nlat: nlat + n2] = tmp[ind2]

#                 count0 = np.full((2*nlat), 0.)
                
#                 if ndim==3:
#                     var_tmp0 = np.full((nh, 2*nlat), -999.)
#                     var_tmp0[:,nlat - n1: nlat] = var_tmp[:,ind1]
#                     var_tmp0[:,nlat: nlat + n2] = var_tmp[:,ind2]
#                 else:
#                     var_tmp0 = np.full((2*nlat), -999.)
#                     var_tmp0[nlat - n1: nlat] = var_tmp[ind1]
#                     var_tmp0[nlat: nlat + n2] = var_tmp[ind2]
#             else:
#                 tmp1 = np.full((2*nlat), -999.)
#                 tmp1[nlat - n1: nlat] = tmp[ind1]
#                 tmp1[nlat: nlat + n2] = tmp[ind2]

#                 count1 = np.full((2*nlat), 0.)
                
#                 if ndim==3:
#                     var_tmp1 = np.full((nh, 2*nlat), -999.)
#                     var_tmp1[:,nlat - n1: nlat] = var_tmp[:,ind1]
#                     var_tmp1[:,nlat: nlat + n2] = var_tmp[:,ind2]
#                 else:
#                     var_tmp1 = np.full((2*nlat), -999.)
#                     var_tmp1[nlat - n1: nlat] = var_tmp[ind1]
#                     var_tmp1[nlat: nlat + n2] = var_tmp[ind2]

#                 tmp0 = np.dstack((tmp0, tmp1))
#                 count0 = np.dstack((count0, count1))
#                 var_tmp0 = np.dstack((var_tmp0, var_tmp1))

#         if ndim == 2:
#             count0[var_tmp0 != -999.] = 1.
            
#         var_tmp0[var_tmp0 == -999.] = np.nan
#         tmp0[tmp0 == -999.] = np.nan
        
#         count00 = nlon*np.nanmean(count0[0], axis = 1)
#         tmp00 = np.nanmean(tmp0[0], axis = 1)
#         if ndim == 3:
#             var_tmp00 = np.nanmean(var_tmp0, axis = 2)
#         else:
#             var_tmp00 = np.nanmean(var_tmp0[0], axis = 1)

#         if j == 0:
#             count_out = count00
#             var_out = var_tmp00
#             ice_out = tmp00
#         else:
#             count_out = np.dstack((count_out, count00))        
#             var_out = np.dstack((var_out, var_tmp00))
#             ice_out = np.dstack((ice_out, tmp00))

#     count_out = count_out[0]
#     ice_out = ice_out[0]
#     if ndim != 3:
#         var_out = var_out[0]
    
#     return var_out, ice_out