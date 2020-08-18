function test_anim(fn)
    scene = Scene()
    record(scene, fn, 1:5) do n
        lines!(scene, 1:5, (1:5) .* n)
    end
end
