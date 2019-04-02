def power_of_two(n):
    """Calculates :math:`2^n`.
    
    Calculates non-negative powers of two via bit shift.
    
    Args:
        n (int): the exponent to raise 2 to.
    Returns:
        int: :math:`2^n`.
    Raises:
       ValueError: if :math:`n \lt 0`.
    """
    
    # left bit shift by n equivalent to 2^n.
    return 1 << n
