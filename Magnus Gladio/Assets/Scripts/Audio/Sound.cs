using UnityEngine.Audio;
using UnityEngine;
using System;

[Serializable]
public class Sound
{

    public string name;
    public AudioClip clip;

    [Range(0, 1)]
    public float volume;
    public float maxVolume;
    [Range(0.1f, 3f)]
    public float pitch;

    [HideInInspector]
    public AudioSource source;

    public bool loop;
}
