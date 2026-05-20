function CodeBlock(el)
  if el.classes:includes("math") then
    return pandoc.Para({pandoc.Math("DisplayMath", el.text)})
  end
end
