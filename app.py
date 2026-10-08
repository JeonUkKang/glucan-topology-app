import streamlit as st
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import io

# =============================================================
# 1. 페이지 설정 및 헤더
# =============================================================
st.set_page_config(page_title="Advanced Glucan Topology Visualizer", layout="wide")

st.markdown("""
    <h2 style='color:#0B3C5D; border-bottom: 2px solid #0B3C5D; padding-bottom: 5px; margin-bottom: 10px;'>
    🧬 Advanced Glucan Topology Visualizer
    </h2>
    <p style='color:#555; font-size: 14px; margin-bottom: 20px;'>
    Data-Driven Macroscopic Topology Rendering Engine
    </p>
""", unsafe_allow_html=True)

# =============================================================
# 2. 프리셋(Preset) 데이터 및 세션 상태(Session State) 동기화 로직
# =============================================================
PRESETS = {
    "[Low DP] Maltotriose": {"mode": "Direct DP Input", "dp": 3, "1,4": 66.7, "t": 33.3},
    "[Low DP] Isomaltotriose": {"mode": "Direct DP Input", "dp": 3, "1,6": 66.7, "t": 33.3},
    "[Low DP] Panose": {"mode": "Direct DP Input", "dp": 3, "1,4": 33.3, "1,6": 33.3, "t": 33.4},
    "[Low DP] Isopanose": {"mode": "Direct DP Input", "dp": 3, "1,4": 33.3, "1,6": 33.3, "t": 33.4},
    "[Low DP] Maltotetraose": {"mode": "Direct DP Input", "dp": 4, "1,4": 75.0, "t": 25.0},
    "[Low DP] DP10-Maltooligosaccharide": {"mode": "Direct DP Input", "dp": 10, "1,4": 90.0, "t": 10.0},
    "[Low DP] DP15-Maltooligosaccharide": {"mode": "Direct DP Input", "dp": 15, "1,4": 93.3, "t": 6.7},
    "[Polymer] High α-(1,2,6) branched Dextran": {"mode": "Calculate from Mw", "mw": 101000.0, "1,6": 34.5, "2,6": 32.4, "t": 33.1},
    "[Polymer] Complex Branched Dextran": {"mode": "Calculate from Mw", "mw": 5470.0, "1,6": 49.8, "1,3": 2.8, "1,4": 0.7, "3,6": 17.6, "2,6": 1.4, "t": 27.7},
    "[Polymer] Alternan": {"mode": "Direct DP Input", "dp": 300, "1,6": 45.1, "1,3": 35.3, "3,6": 9.8, "t": 9.8},
    "[Polymer] Reuteran": {"mode": "Direct DP Input", "dp": 300, "1,4": 45.0, "1,6": 27.0, "4,6": 15.0, "t": 13.0},
    "[Polymer] Mutan": {"mode": "Direct DP Input", "dp": 300, "1,3": 66.7, "1,6": 10.5, "3,6": 11.4, "t": 11.4},
    "[Polymer] Amylopectin": {"mode": "Calculate from Mw", "mw": 100000.0, "1,4": 95.8, "4,6": 4.2, "t": 4.2},
    "[Polymer] Glycogen": {"mode": "Calculate from Mw", "mw": 2270.0, "1,4": 80.7, "4,6": 12.4, "t": 1.1},
    "[Polymer] Pullulan": {"mode": "Calculate from Mw", "mw": 350.0, "1,4": 66.7, "1,6": 33.3}
}

preset_options = [
    "Custom (Manual Input)",
    "-------------------------",
    "[Low DP] Maltotriose",
    "[Low DP] Isomaltotriose",
    "[Low DP] Panose",
    "[Low DP] Isopanose",
    "[Low DP] Maltotetraose",
    "[Low DP] DP10-Maltooligosaccharide",
    "[Low DP] DP15-Maltooligosaccharide",
    "-------------------------",
    "[Polymer] High α-(1,2,6) branched Dextran",
    "[Polymer] Complex Branched Dextran",
    "[Polymer] Alternan",
    "[Polymer] Reuteran",
    "[Polymer] Mutan",
    "[Polymer] Amylopectin",
    "[Polymer] Glycogen",
    "[Polymer] Pullulan"
]

# 초기 변수 세팅
default_keys = {
    "sn": "", "dp_mode": "Direct DP Input", "tgt_dp": 300, "mw_val": 500.0, "disp_glc": 18,
    "t_glc": 0.0, "g12": 0.0, "g13": 0.0, "g14": 0.0, "g16": 0.0, 
    "g26": 0.0, "g36": 0.0, "g46": 0.0, "g236": 0.0
}
for k, v in default_keys.items():
    if k not in st.session_state:
        st.session_state[k] = v

def load_preset():
    sel = st.session_state.preset_selector
    if sel in PRESETS:
        p = PRESETS[sel]
        st.session_state.sn = sel.split("] ")[-1]
        st.session_state.dp_mode = p.get("mode", "Direct DP Input")
        st.session_state.tgt_dp = int(p.get("dp", 300))
        st.session_state.mw_val = float(p.get("mw", 500.0))
        st.session_state.t_glc = float(p.get("t", 0.0))
        st.session_state.g12 = float(p.get("1,2", 0.0))
        st.session_state.g13 = float(p.get("1,3", 0.0))
        st.session_state.g14 = float(p.get("1,4", 0.0))
        st.session_state.g16 = float(p.get("1,6", 0.0))
        st.session_state.g26 = float(p.get("2,6", 0.0))
        st.session_state.g36 = float(p.get("3,6", 0.0))
        st.session_state.g46 = float(p.get("4,6", 0.0))
        st.session_state.g236 = float(p.get("2,3,6", 0.0))
    elif sel == "Custom (Manual Input)":
        st.session_state.sn = ""
        st.session_state.t_glc = 0.0
        st.session_state.g12 = 0.0
        st.session_state.g13 = 0.0
        st.session_state.g14 = 0.0
        st.session_state.g16 = 0.0
        st.session_state.g26 = 0.0
        st.session_state.g36 = 0.0
        st.session_state.g46 = 0.0
        st.session_state.g236 = 0.0

# =============================================================
# 3. UI 레이아웃 및 입력
# =============================================================
st.markdown("#### 1. Sample & Polymer Settings")

st.selectbox("📚 Load Preset:", options=preset_options, key="preset_selector", on_change=load_preset)

col1, col2 = st.columns([2, 1])
with col1:
    st.text_input("Sample Name:", key="sn", placeholder="e.g., Maltotriose, Panose, Mutan")
with col2:
    st.write("") 
    st.write("")
    show_bracket = st.checkbox("Show Brackets (n)", value=True, help="Displays the estimated repeating unit range (n).")

col3, col4, col5, col6 = st.columns(4)
with col3:
    st.selectbox("DP Mode:", options=["Calculate from Mw", "Direct DP Input"], key="dp_mode")
with col4:
    st.number_input("Target DP:", step=10, key="tgt_dp")
with col5:
    st.number_input("Mw (kDa):", step=10.0, key="mw_val")
with col6:
    st.number_input("Display Glc (ea):", min_value=5, max_value=50, step=1, key="disp_glc", help="Will automatically scale down for DP <= 20.")

st.markdown("<hr style='margin: 10px 0px;'>", unsafe_allow_html=True)
st.markdown("#### 2. Linkage Stoichiometry (%)")

left_col, right_col = st.columns(2)
with left_col:
    st.markdown("**[ Linear & Terminal ]**")
    l_c1, l_c2 = st.columns(2)
    with l_c1:
        st.number_input("t-Glc (%):", step=0.1, key="t_glc")
        st.number_input("1,3-Glc (%):", step=0.1, key="g13")
        st.number_input("1,6-Glc (%):", step=0.1, key="g16")
    with l_c2:
        st.number_input("1,2-Glc (%):", step=0.1, key="g12")
        st.number_input("1,4-Glc (%):", step=0.1, key="g14")

with right_col:
    st.markdown("**[ Branching Points ]**")
    r_c1, r_c2 = st.columns(2)
    with r_c1:
        st.number_input("2,6-Glc (%):", step=0.1, key="g26")
        st.number_input("4,6-Glc (%):", step=0.1, key="g46")
    with r_c2:
        st.number_input("3,6-Glc (%):", step=0.1, key="g36")
        st.number_input("2,3,6-Glc (%):", step=0.1, key="g236")

# =============================================================
# 4. 렌더링 엔진 코어 로직
# =============================================================
def draw_bracket(ax, p_start, p_end, normal_vec, length=0.85, lw=1.8, color='black'):
    v_norm = normal_vec / np.linalg.norm(normal_vec)
    v_along = (p_end - p_start) / np.linalg.norm(p_end - p_start)
    tick = v_along * 0.35
    for p in [p_start, p_end]:
        dir_tick = tick if p is p_start else -tick
        pt_a = p + v_norm * length + dir_tick
        pt_b = p + v_norm * length
        pt_c = p - v_norm * length
        pt_d = p - v_norm * length + dir_tick
        ax.plot([pt_a[0], pt_b[0], pt_c[0], pt_d[0]], [pt_a[1], pt_b[1], pt_c[1], pt_d[1]], color=color, lw=lw, zorder=5)

def put_text(ax, pt1, pt2, text):
    mid = (pt1 + pt2) / 2.0
    v = pt2 - pt1
    norm = np.linalg.norm(v)
    if norm == 0: return
    u = v / norm
    n_vec = np.array([u[1], -u[0]])
    if n_vec[1] > 0.05: n_vec = -n_vec
    if abs(n_vec[1]) <= 0.05 and n_vec[0] < 0: n_vec = -n_vec
    offset_dist = 0.32
    pos = mid + n_vec * offset_dist
    ax.text(pos[0], pos[1], text, fontsize=6.5, fontweight='bold', ha='center', va='center')

def extract_topology_data(t_g, g16, g13, g14, g12, g36, g46, g26, g236):
    linear_dict = {"1,4": g14, "1,6": g16, "1,3": g13, "1,2": g12}
    linear_sorted = sorted([(k, v) for k, v in linear_dict.items() if v > 0], key=lambda x: x[1], reverse=True)
    branch_ratio = g36 + g46 + g26 + (g236 * 2)
    if not linear_sorted:
        backbone_str = "1,4-Glc (100.0%)"
        primary_link = "1,4"
    else:
        backbone_str = " + ".join([f"{k}-Glc ({v:.1f}%)" for k, v in linear_sorted])
        primary_link = linear_sorted[0][0]
    return {"primary_link": primary_link, "backbone_str": backbone_str, "db": branch_ratio}

def render_stoichiometric_models(sample_name, rep_dp, show_bracket, num_bb, t_g, g16, g13, g14, g12, g36, g46, g26, g236):
    branch_ratio = g36 + g46 + g26 + (g236 * 2)
    linear_dict = {"1,4": g14, "1,6": g16, "1,3": g13, "1,2": g12}
    branch_dict = {"3,6-Glc": g36, "4,6-Glc": g46, "2,6-Glc": g26, "2,3,6-Glc": g236}
    
    linear_sorted = sorted([(k, v) for k, v in linear_dict.items() if v > 0], key=lambda x: x[1], reverse=True)
    if not linear_sorted: linear_sorted = [("1,4", 100.0)]
    primary_link = linear_sorted[0][0]
    primary_pct = linear_sorted[0][1]

    branch_sorted = sorted([(k, v) for k, v in branch_dict.items() if v > 0], key=lambda x: x[1], reverse=True)
    main_br_name = branch_sorted[0][0] if branch_sorted else "t-Glc"

    if primary_link == "1,3": branch_bond = "1,6"
    elif primary_link == "1,6": branch_bond = "1,3"
    else: branch_bond = "1,6" 

    total_bonds = num_bb - 1
    total_linear_pct = sum([v for k, v in linear_sorted])
    base_len = 5 if primary_pct > 85.0 else (4 if primary_pct > 60.0 else 3)
    avg_chain_len = 100.0 / max(0.1, t_g) if t_g > 0 else 100.0
    use_dashed_extension = avg_chain_len > (base_len + 1) and rep_dp > 20

    alt_links = [primary_link] * total_bonds
    if len(linear_sorted) > 1:
        min_pct = linear_sorted[-1][1]
        is_periodic = True
        pattern_counts = []
        for k, v in linear_sorted:
            ratio = v / min_pct
            if abs(ratio - round(ratio)) > 0.15: 
                is_periodic = False
                break
            pattern_counts.append((k, int(round(ratio))))
        
        if is_periodic:
            pattern = []
            for k, c in pattern_counts:
                pattern.extend([k] * c)
            alt_links = [pattern[i % len(pattern)] for i in range(total_bonds)]
        else:
            counts = [{"link": k, "count": int(np.round(total_bonds * (v / total_linear_pct)))} for k, v in linear_sorted]
            diff = total_bonds - sum(c["count"] for c in counts)
            if diff != 0: counts[0]["count"] += diff

            alt_links = [None] * total_bonds
            idx_list = list(range(total_bonds))
            for c in counts:
                if c["count"] <= 0: continue
                step = len(idx_list) / c["count"]
                for i in range(c["count"]):
                    idx = int(i * step + step / 2)
                    alt_links[idx_list[idx]] = c["link"]
                idx_list = [x for x in range(total_bonds) if alt_links[x] is None]

    branch_freq = branch_ratio / max(1.0, total_linear_pct + branch_ratio)
    num_branches = int(np.round(num_bb * branch_freq))
    if num_branches < 1 and branch_ratio > 0: num_branches = 1

    safe_spots = [i for i in range(1, total_bonds) if alt_links[i-1] == alt_links[i] and alt_links[i] != branch_bond]
    if len(safe_spots) < num_branches:
        safe_spots = [i for i in range(1, total_bonds) if alt_links[i] != branch_bond]

    b_indices = []
    if num_branches > 0 and safe_spots:
        actual_branches = min(num_branches, len(safe_spots))
        step = len(safe_spots) / actual_branches
        for i in range(actual_branches):
            idx = safe_spots[int(i * step + step / 2)] 
            b_indices.append(idx)
    b_indices = sorted(list(set(b_indices)))

    branch_specs = []
    for k_i, b_idx in enumerate(b_indices):
        root_link = branch_bond 
        side_len = min(base_len, max(1, rep_dp - 2)) if rep_dp <= 20 else base_len
        side_links = [primary_link] * (side_len - 1) 
        sub_links = []
        if k_i % 2 == 0 and len(side_links) >= 1 and rep_dp > 10:
            sub_root = branch_bond   
            sub_next = primary_link  
            sub_len = max(1, side_len - 2)
            sub_links = [sub_root] + [sub_next] * sub_len
        branch_specs.append({'root': root_link, 'side': side_links, 'sub_links': sub_links})

    def get_vec(l_type, tier=0):
        if primary_link in ["1,4", "1,3", "1,2"]: 
            if tier == 0: return np.array([1.0, -1.2]) if l_type == "1,6" else np.array([1.4, 0.0]) 
            elif tier == 1: return np.array([1.4, 0.0]) if l_type == primary_link else np.array([1.0, 1.2]) 
            elif tier == 2: return np.array([1.4, 0.0]) if l_type == primary_link else np.array([1.0, 1.2]) 
        else: 
            if tier == 0: return np.array([1.0, 1.0]) if l_type == "1,6" else np.array([1.4, 0.0])
            elif tier == 1: return np.array([0.0, 1.4]) if l_type == "1,6" else np.array([-1.0, 1.0]) 
            elif tier == 2: return np.array([-1.0, 1.0]) if l_type == "1,6" else np.array([-1.4, 0.0]) 
        return np.array([1.4, 0.0])

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(22, 11.5))
    fig.patch.set_facecolor('white')

    models = [("Topology I: Single-Tier Branching\n(Comb-like)", ax1, 1), ("Topology II: Multi-Tier Branching\n(Hierarchical / Ramified)", ax2, 2)]

    for title, ax, m_type in models:
        ax.set_title(title, fontsize=14, fontweight='bold', pad=15, loc='center', color='#0B3C5D')
        ax.set_aspect('equal')
        ax.axis('off')
        
        all_coords = []
        t_glc_count = 0 

        bb_coords = [np.array([0.0, 0.0])]
        for i in range(total_bonds):
            bb_coords.append(bb_coords[-1] + get_vec(alt_links[i], tier=0))
        bb_coords = np.array(bb_coords)
        all_coords.extend(bb_coords)

        for i in range(total_bonds):
            ax.plot([bb_coords[i,0], bb_coords[i+1,0]], [bb_coords[i,1], bb_coords[i+1,1]], color='black', lw=2.2, zorder=1)
            put_text(ax, bb_coords[i], bb_coords[i+1], alt_links[i])

        if rep_dp > 20:
            start_v, end_v = get_vec(alt_links[0], tier=0) * 0.8, get_vec(alt_links[-1], tier=0) * 0.8
            ax.plot([bb_coords[0,0] - start_v[0], bb_coords[0,0]], [bb_coords[0,1] - start_v[1], bb_coords[0,1]], color='#333333', lw=2.2, linestyle='--', dashes=(2.0, 1.2), zorder=1)
            ax.plot([bb_coords[-1,0], bb_coords[-1,0] + end_v[0]], [bb_coords[-1,1], bb_coords[-1,1] + end_v[1]], color='#333333', lw=2.2, linestyle='--', dashes=(2.0, 1.2), zorder=1)
            all_coords.extend([bb_coords[0] - start_v, bb_coords[-1] + end_v])

        bb_name = f"{primary_link}-Glc"
        for i, (bx, by) in enumerate(bb_coords):
            ax.add_patch(plt.Circle((bx, by), 0.32, facecolor='#1f77b4', edgecolor='black', lw=1.8, zorder=3))
            if i in [0, num_bb - 1] and rep_dp > 20: 
                ax.text(bx, by - 0.45, f"{bb_name}\n{primary_pct:.1f}%", fontsize=7.0, ha='center', va='top')
            elif i in b_indices: 
                ax.text(bx, by - 0.45, f"{main_br_name}\n{branch_ratio:.1f}%", fontsize=6.5, ha='center', va='top')

        for k_i, b_idx in enumerate(b_indices):
            px, py = bb_coords[b_idx]
            spec = branch_specs[k_i]
            
            if m_type == 1:
                t1_start = np.array([px, py]) + get_vec(spec['root'], tier=1)
                ax.plot([px, t1_start[0]], [py, t1_start[1]], color='black', lw=2.0, zorder=1)
                ax.add_patch(plt.Circle(t1_start, 0.32, facecolor='#1f77b4', edgecolor='black', lw=1.8, zorder=3))
                put_text(ax, np.array([px, py]), t1_start, spec['root'])
                if t_glc_count < 2 and t_g > 0:
                    ax.text(t1_start[0], t1_start[1] + 0.45, f"t-Glc\n{t_g:.1f}%", fontsize=6.5, ha='center', va='bottom')
                    t_glc_count += 1
                all_coords.append(t1_start)
                
            elif m_type == 2:
                t1_start = np.array([px, py]) + get_vec(spec['root'], tier=1)
                ax.plot([px, t1_start[0]], [py, t1_start[1]], color='black', lw=2.0, zorder=1)
                ax.add_patch(plt.Circle(t1_start, 0.32, facecolor='#1f77b4', edgecolor='black', lw=1.8, zorder=3))
                put_text(ax, np.array([px, py]), t1_start, spec['root'])
                
                t1_coords = [t1_start.copy()]
                cur1 = t1_start.copy()
                for l_type in spec['side']:
                    cur1 += get_vec(l_type, tier=1)
                    t1_coords.append(cur1.copy())
                for k in range(len(t1_coords)-1):
                    ax.plot([t1_coords[k][0], t1_coords[k+1][0]], [t1_coords[k][1], t1_coords[k+1][1]], color='black', lw=2.0, zorder=1)
                    ax.add_patch(plt.Circle(t1_coords[k+1], 0.32, facecolor='#1f77b4', edgecolor='black', lw=1.8, zorder=3))
                    put_text(ax, t1_coords[k], t1_coords[k+1], spec['side'][k])
                all_coords.extend(t1_coords)
                t1_text_pos = t1_coords[-1]

                t2_text_pos = None
                if spec['sub_links'] and rep_dp > 10:
                    delay_idx = min(2, len(t1_coords) - 1)
                    t2_root_start = t1_coords[delay_idx] 
                    t2_start = t2_root_start + get_vec(spec['sub_links'][0], tier=2)
                    ax.plot([t2_root_start[0], t2_start[0]], [t2_root_start[1], t2_start[1]], color='black', lw=2.0, zorder=1)
                    ax.add_patch(plt.Circle(t2_start, 0.32, facecolor='#1f77b4', edgecolor='black', lw=1.8, zorder=3))
                    put_text(ax, t2_root_start, t2_start, spec['sub_links'][0])
                    ax.text(t2_root_start[0] + 0.3, t2_root_start[1] - 0.2, f"{main_br_name}\n{branch_ratio:.1f}%", fontsize=6.0, ha='left', va='top')
                    
                    t2_coords = [t2_start.copy()]
                    cur2 = t2_start.copy()
                    for l_type in spec['sub_links'][1:]:
                        cur2 += get_vec(l_type, tier=2)
                        t2_coords.append(cur2.copy())
                    for k in range(len(t2_coords)-1):
                        ax.plot([t2_coords[k][0], t2_coords[k+1][0]], [t2_coords[k][1], t2_coords[k+1][1]], color='black', lw=2.0, zorder=1)
                        ax.add_patch(plt.Circle(t2_coords[k+1], 0.32, facecolor='#1f77b4', edgecolor='black', lw=1.8, zorder=3))
                        put_text(ax, t2_coords[k], t2_coords[k+1], spec['sub_links'][1:][k])
                    all_coords.extend(t2_coords)
                    t2_text_pos = t2_coords[-1]

                if use_dashed_extension:
                    ex_vec1 = get_vec(primary_link, tier=1)
                    d1_1, d1_2 = t1_coords[-1] + ex_vec1 * 0.4, t1_coords[-1] + ex_vec1 * 1.5
                    ax.plot([t1_coords[-1][0], d1_1[0]], [t1_coords[-1][1], d1_1[1]], color='black', lw=2.0, zorder=1)
                    ax.plot([d1_1[0], d1_2[0]], [d1_1[1], d1_2[1]], color='#333333', lw=2.0, linestyle='--', dashes=(2.0, 1.2), zorder=1)
                    t1_text_pos = d1_2

                    if t2_text_pos is not None:
                        ex_vec2 = get_vec(primary_link, tier=2)
                        d2_1, d2_2 = t2_coords[-1] + ex_vec2 * 0.4, t2_coords[-1] + ex_vec2 * 1.5
                        ax.plot([t2_coords[-1][0], d2_1[0]], [t2_coords[-1][1], d2_1[1]], color='black', lw=2.0, zorder=1)
                        ax.plot([d2_1[0], d2_2[0]], [d2_1[1], d2_2[1]], color='#333333', lw=2.0, linestyle='--', dashes=(2.0, 1.2), zorder=1)
                        t2_text_pos = d2_2

                if t_g > 0 and t_glc_count < 2:
                    ax.text(t1_text_pos[0], t1_text_pos[1] + 0.35, f"t-Glc\n{t_g:.1f}%", fontsize=6.5, ha='center', va='bottom')
                    t_glc_count += 1
                    if t2_text_pos is not None:
                        ax.text(t2_text_pos[0], t2_text_pos[1] + 0.35, f"t-Glc\n{t_g:.1f}%", fontsize=6.5, ha='center', va='bottom')
                        t_glc_count += 1

        all_x, all_y = [pt[0] for pt in all_coords], [pt[1] for pt in all_coords]
        ax.set_xlim(min(all_x) - 2.5, max(all_x) + 3.0)
        ax.set_ylim(min(all_y) - 3.5, max(all_y) + 4.5)

        if show_bracket and len(bb_coords) >= 3 and rep_dp > 5:
            p_start, p_end = (bb_coords[1] + bb_coords[2]) / 2.0, (bb_coords[-2] + bb_coords[-3]) / 2.0
            n_vec_b = np.array([0, -1]) if primary_link in ["1,4", "1,3", "1,2"] else np.array([1, -1]) 
            draw_bracket(ax, p_start, p_end, n_vec_b, length=0.85, lw=1.8, color='black')
            label_pos = p_end + (n_vec_b / np.linalg.norm(n_vec_b)) * 1.5 + np.array([0.5, 0.0])
            rep_n_low = max(1, int(np.round((rep_dp * (primary_pct / 100.0)) / num_bb)))
            rep_n_high = max(1, int(np.round((rep_dp * (primary_pct / 100.0)) / (num_bb * 0.65))))
            ax.text(label_pos[0], label_pos[1], f"$n \\approx {rep_n_low}-{rep_n_high}$", fontsize=9, fontweight='bold', ha='left', va='center')

    plt.suptitle(f"Topological Spectrum for {sample_name if sample_name else 'Oligosaccharide'} (DP = {rep_dp})", fontsize=16, fontweight='bold', y=0.98, color='#0B3C5D')
    plt.tight_layout(rect=[0, 0, 1, 0.92])
    
    return fig

# =============================================================
# 5. 버튼 이벤트 및 출력
# =============================================================
st.markdown("<br>", unsafe_allow_html=True)
if st.button("🚀 Generate Theoretical Topologies", use_container_width=True):
    
    # 세션 상태(Session State)에서 값 불러오기
    rep_dp = int(np.round((st.session_state.mw_val * 1000.0) / 162.14)) if st.session_state.dp_mode == "Calculate from Mw" else int(st.session_state.tgt_dp)
    if rep_dp < 1: rep_dp = 3

    # 🚨 저분자(DP 3~20) 자동 적응 로직
    if rep_dp <= 20:
        num_bb = rep_dp
    else:
        num_bb = int(st.session_state.disp_glc) if st.session_state.disp_glc > 0 else 18
    
    t_g, g12, g13, g14, g16 = st.session_state.t_glc, st.session_state.g12, st.session_state.g13, st.session_state.g14, st.session_state.g16
    g26, g36, g46, g236 = st.session_state.g26, st.session_state.g36, st.session_state.g46, st.session_state.g236
    s_name = st.session_state.sn

    total_sum = t_g + g12 + g13 + g14 + g16 + g26 + g36 + g46 + g236
    eval_data = extract_topology_data(t_g, g16, g13, g14, g12, g36, g46, g26, g236)
    
    avg_cl = f"{100.0 / t_g:.1f}" if t_g > 0 else "N/A"
    branch_dict = {"3,6-Glc": g36, "4,6-Glc": g46, "2,6-Glc": g26, "2,3,6-Glc": g236}
    branch_sorted = sorted([(k, v) for k, v in branch_dict.items() if v > 0], key=lambda x: x[1], reverse=True)
    main_br_str = " + ".join([f"{k} ({v:.1f}%)" for k, v in branch_sorted]) if branch_sorted else "None"

    # Report 컨테이너
    with st.container():
        st.markdown("---")
        st.markdown("### 📋 Topological Structural Report")
        colA, colB = st.columns(2)
        with colA:
            st.write(f"**Sample Name:** {s_name if s_name else 'Oligosaccharide / Polymer'}")
            st.write(f"**Target DP (Mw):** {rep_dp} ({st.session_state.mw_val} kDa)")
            st.write(f"**Stoichiometry Sum:** {total_sum:.1f}%")
            st.write(f"**Linear Backbone(s):** {eval_data['backbone_str']}")
        with colB:
            st.write(f"**Main Branch Linkage:** {main_br_str}")
            st.write(f"**Degree of Branching (DB):** {eval_data['db']:.1f}%")
            st.write(f"**Terminal Units:** {t_g:.1f}%")
            st.write(f"**Average Chain Length (CL):** {avg_cl}")
        
        st.info("▶ The Data-Driven Engine visualizes theoretically possible structural candidates.\n\n▶ Researchers should independently select the final topological model (Topology I or II) based on DB, steric hindrance, and structural assays.")
        st.markdown("---")

    # 렌더링 및 출력
    fig = render_stoichiometric_models(s_name, rep_dp, show_bracket, num_bb, t_g, g16, g13, g14, g12, g36, g46, g26, g236)
    st.pyplot(fig)

    # 고해상도 파일 다운로드
    png_buffer = io.BytesIO()
    fig.savefig(png_buffer, format='png', dpi=600, bbox_inches='tight')
    png_buffer.seek(0)
    
    pdf_buffer = io.BytesIO()
    fig.savefig(pdf_buffer, format='pdf', bbox_inches='tight')
    pdf_buffer.seek(0)
    
    st.markdown("<br>", unsafe_allow_html=True)
    down_col1, down_col2, _ = st.columns([1, 1, 2])
    with down_col1:
        st.download_button(
            label="📥 Download as PNG (300dpi)",
            data=png_buffer,
            file_name=f"{s_name.replace(' ', '_') if s_name else 'topology'}.png",
            mime="image/png",
            use_container_width=True
        )
    with down_col2:
        st.download_button(
            label="📥 Download as PDF (Vector)",
            data=pdf_buffer,
            file_name=f"{s_name.replace(' ', '_') if s_name else 'topology'}.pdf",
            mime="application/pdf",
            use_container_width=True
        )
    plt.close(fig)
