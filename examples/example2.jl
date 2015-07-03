using Pathogen

create_tree([generate_seq(200, 0.25, 0.25, 0.25, 0.25),
             generate_seq(200, 0.25, 0.25, 0.25, 0.25),
             generate_seq(200, 0.25, 0.25, 0.25, 0.25),
             generate_seq(200, 0.25, 0.25, 0.25, 0.25)],
            [10., 9., 8., 7.])
