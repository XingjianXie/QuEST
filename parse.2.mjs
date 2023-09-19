// read xml file from all.xml in js
import fs from 'fs';
import { XMLParser } from 'fast-xml-parser';
import { exec } from 'child_process';
import { promisify } from 'util';

function defmapper(def) {
    return def.map(member => {
        if (member.references && !Array.isArray(member.references)) {
            member.references = [member.references];
        }
        return {
            name: member.name,
            references: member.references ? member.references.filter(
                ref =>
                ref['@_refid'].startsWith('QuEST__common_8c') ||
                ref['@_refid'].startsWith('QuEST__internal_8h')
                ).map(
                    ref => {
                        return {
                            name: ref['#text'],
                            type: ref['@_refid'].startsWith('QuEST__common_8c') ? 'common' : 'internal'
                        }
                    }
                ) : []
        }
    })
}

function addref(common, gf, fr, set) {
    for (const ref of fr) {
        if (ref.type === 'internal') {
            gf.push({name: ref.name, type: 'internal'});
        } else {
            const elem = {name: ref.name, type: 'common', dependencies: []};
            if (set.has(ref.name)) {
                elem.dependencies = 'recursive'
                gf.push(elem);
            } else {
                set.add(ref.name);
                addref(common, elem.dependencies, common.find(c => c.name === ref.name)?.references ?? [], set)
                gf.push(elem);
                set.delete(ref.name);
            }
        }
    }
}

async function main() {
    const xml = fs.readFileSync('./doc/xml/all.xml', 'utf-8');
    const parser = new XMLParser({
        ignoreAttributes: false,
    });
    const obj = parser.parse(xml);
    const compounddef = obj.doxygen.compounddef.reduce((acc, cur) => {
        acc[cur.compoundname] = cur
        return acc
    }, {});


    const questh = compounddef['QuEST.h'].sectiondef.find(def => def['@_kind'] === 'func').member;

    const cmd = `nm -gU build/QuEST/libQuEST.dylib | c++filt | grep ' T ' | awk '{print $3}'`

    // execute command using es module
    const { stdout, stderr } = await promisify(exec)(cmd);

    const pub = stdout.trim().split('\n');

    for (const def of questh) {
        if (pub.find(p => p === "_" + def.name)) {
            continue;
        } else {
            console.log("Missing public function: " + def.name)
        }
    }
}

main();
