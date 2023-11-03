from __future__ import annotations

import numpy as np
import numpy.typing as npt
from pymatching import Matching
from scipy.sparse import block_diag, csr_matrix, eye, hstack
from utils.simulation_utils import (
    get_virtual_check_init_vals,
)


def build_multiround_pcm(pcm, repetitions, format="csr") -> npt.NDArray[int]:
    """Builds the multiround parity-check matrix as described in the paper.

    Each row corresponds to a round of measurements, the matrix for r repetitions has the form
    H3D := (H_diag | id_r), where H_diag is the r-diagonal matrix containing the parity-check matrix H and
    id_diag is a staircase block matrix containing two m-identity matrices in each row and column except the first one.
    """
    if not isinstance(pcm, csr_matrix):
        pcm = csr_matrix(pcm)

    pcm_rows, pcm_cols = pcm.shape

    # Construct the block of PCMs
    H_3DPCM = block_diag([pcm] * (repetitions + 1), format=format)

    # Construct the block of identity matrices
    H_3DID_diag = block_diag([eye(pcm_rows, format=format)] * (repetitions + 1), format=format)

    # Construct the block of identity matrices
    H_3DID_offdiag = eye(pcm_rows * (repetitions + 1), k=-pcm_rows, format=format)

    # Construct the block of identity matrices
    H_3DID = H_3DID_diag + H_3DID_offdiag

    # # hstack the two blocks
    return hstack([H_3DPCM, H_3DID], format=format)


def move_syndrome(syndrome, data_type=np.int32) -> npt.NDArray[int]:
    """Slides the window one region up, i.e., the syndrome of the first half is overwritten by the second half."""
    T = int(syndrome.shape[1] / 2)  # number of rounds in each region

    # zero likes syndrome
    new_syndrome = np.zeros(syndrome.shape, dtype=data_type)

    # move the syndromes from the Tth round to the 0th round
    new_syndrome[:, :T] = syndrome[:, T:]
    return new_syndrome


def get_updated_decoder(decoding_method: str, decoder, new_channel, H3D=None) -> object:
    """Updates the decoder with the new channel information and returns the updated decoder object."""
    if decoding_method == "bposd":
        decoder.update_channel_probs(new_channel)
        return decoder
    elif decoding_method == "matching":
        weights = np.clip(
            np.log((1 - new_channel) / new_channel),
            a_min=-16777215,
            a_max=16777215,
        )
        return Matching(H3D, weights=weights)
    else:
        msg = "Unknown decoding method"
        raise ValueError(msg, decoding_method)


def decode_multiround(
    syndrome: npt.NDArray[int],
    H: npt.NDArray[int],
    decoder,
    channel_probs: npt.NDArray[int],  # needed for matching decoder does not have an update weights method
    repetitions: int,
    last_round=False,
    analog_syndr=None,
    check_block_size: int = 0,
    sigma: float = 0.0,
    H3D: npt.NDArray[int] = None,  # needed for matching decoder
    decoding_method: str = "bposd",  # bposd or matching
) -> tuple[npt.NDArray[int], npt.NDArray[int], npt.NDArray[int], int]:
    """Overlapping window decoding.
    First, we compute the difference syndrome from the recorded syndrome of each measurement round for all measurement
    rounds of the current window (consisting of two regions with equal size).
    Then, we apply the correction returned from the decoder on the first region (commit region).
    We then propagate the syndrome through the whole window (i.e., to the end of the second region).
    """
    analog_tg = analog_syndr is not None
    # convert syndrome to difference syndrome
    diff_syndrome = syndrome.copy()
    diff_syndrome[:, 1:] = (syndrome[:, 1:] - syndrome[:, :-1]) % 2
    bp_iter = 0

    region_size = repetitions // 2  # assumes repetitions is even

    if analog_tg:
        # If we have analog information, we use it to initialize the time-like syndrome nodes, which are defined
        # in the block of the H3D matrix after the diagonal H block.
        analog_init_vals = get_virtual_check_init_vals(analog_syndr.flatten("F"), sigma)

        new_channel = np.hstack((channel_probs[:check_block_size], analog_init_vals))

        # in the last round, we have a perfect syndrome round to make sure we're in the codespace
        if last_round:
            new_channel[-H.shape[0] :] = 1e-15

        decoder = get_updated_decoder(decoding_method, decoder, new_channel, H3D)

    else:
        if last_round:
            new_channel = np.copy(channel_probs)
            new_channel[-H.shape[0] :] = 1e-15

            decoder = get_updated_decoder(decoding_method, decoder, new_channel, H3D)

    decoded = decoder.decode(diff_syndrome.flatten("F"))

    if decoding_method == "bposd":
        bp_iter = decoder.iter

    # extract space correction, first repetitions * n entires
    space_correction = decoded[: H.shape[1] * repetitions].reshape((repetitions, H.shape[1])).T
    # extract time correction

    if last_round is False:
        # this corresponds to the decoding on the second block of the H3D matrix
        time_correction = decoded[H.shape[1] * repetitions :].reshape((repetitions, H.shape[0])).T

        # append time correction with zeros
        time_correction = np.hstack((time_correction, np.zeros((H.shape[0], 1), dtype=np.int32)))

        # correct only in the commit region
        decoded = (np.cumsum(space_correction, 1) % 2)[:, region_size - 1]

        #  get the syndrome according to the correction
        corr_syndrome = (H @ decoded) % 2

        # propagate the syndrome correction through the tentative region
        syndrome[:, region_size:] = ((syndrome[:, region_size:] + corr_syndrome[:, None]) % 2).astype(np.int32)

        # apply the time correction of round region_size - 1 to the syndrome at the beginning of the tentative region
        syndrome[:, region_size] = ((syndrome[:, region_size] + time_correction[:, region_size - 1]) % 2).astype(
            np.int32
        )

    else:
        # correct in the commit and tentative region as the last round stabilizer is perfect
        decoded = (np.cumsum(space_correction, 1) % 2)[:, -1]

    return decoded.astype(np.int32), syndrome, analog_syndr, bp_iter
