--------------------------------------------------------------------------------
-- File         : axis_pixel_extractor.vhd
-- Project      : PALEOCCENE FPGA preprocessing pipeline
-- Description  : Tags each incoming AXI-Stream pixel with its (x, y)
--                coordinate within the current frame. Generates a frame_sync
--                pulse aligned with the first pixel of every frame.
--
--                The AXI-Stream contract assumed here:
--                  - tdata  : 16-bit grayscale pixel
--                  - tvalid : pixel valid
--                  - tuser  : asserted on the first pixel of a frame (SOF)
--                  - tlast  : asserted on the last pixel of a line (EOL)
--                  - tready : always '1' (downstream is pixel-synchronous)
--
--                tuser takes priority over tlast: if both are asserted on
--                the same cycle, that pixel is treated as the start of a
--                new frame and coordinates are reset to (0, 0).
--
-- Notes        : All outputs are registered in a single synchronous process,
--                so pixel_val_out, x_coord, y_coord, data_valid, and
--                frame_sync are always in lockstep.
--
--                The tuser coordinate bug (first pixel of a new frame
--                getting stale coordinates from the previous frame) is
--                fixed by using variables for x/y so that the reset on
--                tuser is evaluated immediately within the process.
--------------------------------------------------------------------------------

library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;

entity axis_pixel_extractor is
    generic (
        FRAME_WIDTH  : positive := 4096;   -- ORCA-Quest 2 full width
        FRAME_HEIGHT : positive := 2304;   -- ORCA-Quest 2 full height
        PIXEL_WIDTH  : positive := 16;     -- grayscale bit depth
        COORD_WIDTH  : positive := 12      -- enough for 4096 (2^12)
    );
    port (
        clk           : in  std_logic;
        rstn          : in  std_logic;

        -- AXI-Stream input
        s_axis_tdata  : in  std_logic_vector(PIXEL_WIDTH-1 downto 0);
        s_axis_tvalid : in  std_logic;
        s_axis_tuser  : in  std_logic;     -- start of frame
        s_axis_tlast  : in  std_logic;     -- end of line
        s_axis_tready : out std_logic;     -- always '1' here

        -- Tagged pixel output (registered, all signals in lockstep)
        pixel_val_out : out std_logic_vector(PIXEL_WIDTH-1 downto 0);
        x_coord       : out std_logic_vector(COORD_WIDTH-1 downto 0);
        y_coord       : out std_logic_vector(COORD_WIDTH-1 downto 0);
        data_valid    : out std_logic;
        frame_sync    : out std_logic      -- pulse on first pixel of frame
    );
end entity axis_pixel_extractor;

architecture rtl of axis_pixel_extractor is

    -- Internal next-position counters (registered between clocks)
    signal next_x_r : unsigned(COORD_WIDTH-1 downto 0) := (others => '0');
    signal next_y_r : unsigned(COORD_WIDTH-1 downto 0) := (others => '0');

begin

    -- Always ready to accept pixels (pure pixel-synchronous downstream)
    s_axis_tready <= '1';

    --------------------------------------------------------------------------
    -- Main coordinate-tagging process
    --
    -- Variables x_v and y_v hold the coordinate of the *current* pixel being
    -- tagged this cycle. They start from the registered next_x_r/next_y_r,
    -- but are overridden immediately to (0,0) if tuser is asserted. This
    -- avoids the stale-coordinate bug where the first pixel of a new frame
    -- would otherwise be tagged with the previous frame's last position.
    --------------------------------------------------------------------------
    process(clk)
        variable x_v : unsigned(COORD_WIDTH-1 downto 0);
        variable y_v : unsigned(COORD_WIDTH-1 downto 0);
    begin
        if rising_edge(clk) then
            if rstn = '0' then
                -- Synchronous reset
                next_x_r      <= (others => '0');
                next_y_r      <= (others => '0');
                pixel_val_out <= (others => '0');
                x_coord       <= (others => '0');
                y_coord       <= (others => '0');
                data_valid    <= '0';
                frame_sync    <= '0';
            else
                -- Default: no valid pixel this cycle
                data_valid <= '0';
                frame_sync <= '0';

                if s_axis_tvalid = '1' then
                    -- Determine coordinates of THIS pixel.
                    -- tuser priority over tlast (matches AXI-Stream spec).
                    if s_axis_tuser = '1' then
                        x_v := (others => '0');
                        y_v := (others => '0');
                        frame_sync <= '1';
                    else
                        x_v := next_x_r;
                        y_v := next_y_r;
                    end if;

                    -- Register outputs for this pixel
                    pixel_val_out <= s_axis_tdata;
                    x_coord       <= std_logic_vector(x_v);
                    y_coord       <= std_logic_vector(y_v);
                    data_valid    <= '1';

                    -- Compute coordinates of the NEXT pixel
                    if s_axis_tlast = '1' then
                        -- End of line: wrap x, advance y (with clamp)
                        x_v := (others => '0');
                        if y_v < to_unsigned(FRAME_HEIGHT-1, COORD_WIDTH) then
                            y_v := y_v + 1;
                        else
                            -- Safety clamp: malformed stream missing tuser.
                            -- Hold at last row; tuser will reset properly.
                            y_v := to_unsigned(FRAME_HEIGHT-1, COORD_WIDTH);
                        end if;
                    else
                        -- Mid-line: advance x (with clamp)
                        if x_v < to_unsigned(FRAME_WIDTH-1, COORD_WIDTH) then
                            x_v := x_v + 1;
                        else
                            -- Safety clamp: malformed stream missing tlast.
                            x_v := to_unsigned(FRAME_WIDTH-1, COORD_WIDTH);
                        end if;
                    end if;

                    -- Store for next cycle
                    next_x_r <= x_v;
                    next_y_r <= y_v;
                end if;
                -- If tvalid = '0', counters and outputs are held (registered)
            end if;
        end if;
    end process;

end architecture rtl;