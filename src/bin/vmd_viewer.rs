use macroquad::prelude::*;
use sang_md::visualization::{load_trajectory, TrajectoryData};

const DEFAULT_FPS: f32 = 24.0;

fn window_conf() -> Conf {
    Conf {
        window_title: "FerrumMD VMD-style Viewer".to_string(),
        window_width: 1280,
        window_height: 720,
        ..Default::default()
    }
}

fn parse_args() -> Result<(String, Option<String>), String> {
    let mut args = std::env::args().skip(1);
    let gro_path = args.next().ok_or_else(|| {
        "usage: cargo run --features viewer --bin vmd_viewer -- <structure.gro> [trajectory.xtc]"
            .to_string()
    })?;
    let xtc_path = args.next();
    Ok((gro_path, xtc_path))
}

fn bounds_for_frame(
    points: &[nalgebra::Vector3<f32>],
) -> (nalgebra::Vector3<f32>, nalgebra::Vector3<f32>) {
    let mut min_v = points[0];
    let mut max_v = points[0];

    for point in points.iter().skip(1) {
        min_v.x = min_v.x.min(point.x);
        min_v.y = min_v.y.min(point.y);
        min_v.z = min_v.z.min(point.z);

        max_v.x = max_v.x.max(point.x);
        max_v.y = max_v.y.max(point.y);
        max_v.z = max_v.z.max(point.z);
    }

    (min_v, max_v)
}

fn compute_scale(data: &TrajectoryData) -> f32 {
    let first_frame = &data.frames[0].positions_nm;
    let (min_v, max_v) = bounds_for_frame(first_frame);
    let extent = max_v - min_v;
    let max_extent = extent.x.max(extent.y).max(extent.z);
    if max_extent <= f32::EPSILON {
        1.0
    } else {
        12.0 / max_extent
    }
}

fn draw_box(box_dims_nm: nalgebra::Vector3<f32>, scale: f32, color: Color) {
    let x = box_dims_nm.x * scale;
    let y = box_dims_nm.y * scale;
    let z = box_dims_nm.z * scale;

    let p000 = vec3(0.0, 0.0, 0.0);
    let p100 = vec3(x, 0.0, 0.0);
    let p010 = vec3(0.0, y, 0.0);
    let p110 = vec3(x, y, 0.0);

    let p001 = vec3(0.0, 0.0, z);
    let p101 = vec3(x, 0.0, z);
    let p011 = vec3(0.0, y, z);
    let p111 = vec3(x, y, z);

    for (a, b) in [
        (p000, p100),
        (p000, p010),
        (p000, p001),
        (p100, p110),
        (p100, p101),
        (p010, p110),
        (p010, p011),
        (p110, p111),
        (p001, p101),
        (p001, p011),
        (p101, p111),
        (p011, p111),
    ] {
        draw_line_3d(a, b, color);
    }
}

#[macroquad::main(window_conf)]
async fn main() {
    let (gro_path, xtc_path) = match parse_args() {
        Ok(values) => values,
        Err(message) => {
            eprintln!("{message}");
            return;
        }
    };

    let data = match load_trajectory(&gro_path, xtc_path.as_deref()) {
        Ok(data) => data,
        Err(message) => {
            eprintln!("trajectory load failed: {message}");
            return;
        }
    };

    let mut frame_idx = 0usize;
    let mut playing = true;
    let mut fps = DEFAULT_FPS;
    let mut elapsed = 0.0f32;
    let mut yaw = -0.6f32;
    let mut pitch = 0.6f32;
    let mut radius = 24.0f32;

    let scale = compute_scale(&data);

    loop {
        let dt = get_frame_time();
        clear_background(BLACK);

        if is_key_pressed(KeyCode::Space) {
            playing = !playing;
        }
        if is_key_pressed(KeyCode::Right) {
            frame_idx = (frame_idx + 1) % data.frames.len();
            playing = false;
        }
        if is_key_pressed(KeyCode::Left) {
            frame_idx = frame_idx.checked_sub(1).unwrap_or(data.frames.len() - 1);
            playing = false;
        }
        if is_key_pressed(KeyCode::Up) {
            fps = (fps + 6.0).min(120.0);
        }
        if is_key_pressed(KeyCode::Down) {
            fps = (fps - 6.0).max(1.0);
        }

        let (mouse_dx, mouse_dy) = mouse_delta_position();
        if is_mouse_button_down(MouseButton::Left) {
            yaw -= mouse_dx * 0.01;
            pitch = (pitch + mouse_dy * 0.01).clamp(-1.4, 1.4);
        }

        let (_, wheel_y) = mouse_wheel();
        radius = (radius - wheel_y * 1.5).clamp(4.0, 100.0);

        if playing && data.frames.len() > 1 {
            elapsed += dt;
            let frame_period = 1.0 / fps;
            while elapsed >= frame_period {
                elapsed -= frame_period;
                frame_idx = (frame_idx + 1) % data.frames.len();
            }
        }

        let frame = &data.frames[frame_idx];
        let mut center = vec3(0.0, 0.0, 0.0);
        for point in &frame.positions_nm {
            center += vec3(point.x * scale, point.y * scale, point.z * scale);
        }
        center /= frame.positions_nm.len() as f32;

        let camera_pos = vec3(
            center.x + radius * yaw.cos() * pitch.cos(),
            center.y + radius * pitch.sin(),
            center.z + radius * yaw.sin() * pitch.cos(),
        );

        set_camera(&Camera3D {
            position: camera_pos,
            up: vec3(0.0, 1.0, 0.0),
            target: center,
            ..Default::default()
        });

        for point in &frame.positions_nm {
            let position = vec3(point.x * scale, point.y * scale, point.z * scale);
            draw_sphere(position, 0.16, None, SKYBLUE);
        }

        if let Some(box_dims) = data.box_dims_nm {
            draw_box(box_dims, scale, GRAY);
        }

        set_default_camera();
        draw_text(
            &format!(
                "FerrumMD Viewer | atoms: {} | frame: {}/{} | step: {} | t = {:.3} ps",
                data.atom_count,
                frame_idx + 1,
                data.frames.len(),
                frame.step,
                frame.time_ps,
            ),
            20.0,
            30.0,
            28.0,
            WHITE,
        );
        draw_text(
            &format!(
                "Space: play/pause | ←/→: step | ↑/↓: FPS ({fps:.0}) | Drag: orbit | Wheel: zoom"
            ),
            20.0,
            58.0,
            24.0,
            LIGHTGRAY,
        );

        next_frame().await;
    }
}
