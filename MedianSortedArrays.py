class Solution:
    def findMedianSortedArrays(self, nums1: List[int], nums2: List[int]) -> float:
        array = sorted(nums1 + nums2)
        
        if len(array) % 2 == 0:
            mid_left = array[int((len(array) / 2) - 1)]
            mid_right = array[int((len(array) / 2))]
            return (mid_left + mid_right) / 2.0
        else:
            mid = array[int(len(array) / 2 - 0.5)]
            return mid
